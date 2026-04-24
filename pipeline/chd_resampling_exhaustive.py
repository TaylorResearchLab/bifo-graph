#!/usr/bin/env python3
"""
chd_resampling_exhaustive.py

Exhaustive in-memory resampling over all C(15,10) = 3,003 possible 10/5 splits
of the 15-gene CHD pool. No per-split files are written. Both PPR operators are
built once; all 3,003 seed vectors are then evaluated via parallel worker
processes. A compact per-split CSV and an aggregated summary JSON are written
on completion.

Usage
-----
    python chd_resampling_exhaustive.py \\
        --kept-edges   test_output/results_kept_edges.csv \\
        --edges-merged test_output/edges_merged.csv \\
        --node-index   test_output/results_node_index.json \\
        --bifo-scores  test_output/pathway_scores_full.csv \\
        --chd-pathways chd_pathway_reference.txt \\
        --out-csv      chd_resampling_results.csv \\
        --out-json     chd_resampling_summary.json \\
        --n-cores      8

Runtime
-------
    ~90 ms/split per core.  Scales near-linearly with --n-cores.
    Approximate wall-clock times for all 3,003 splits:
        1 core  →  ~4.5 min
        4 cores →  ~70 sec
        8 cores →  ~35 sec
       16 cores →  ~18 sec

Parallelism
-----------
    Workers are spawned via ProcessPoolExecutor with an initializer that
    deserialises the sparse PPR operators (passed as CSR component arrays)
    into each worker's address space. Chunks of splits are distributed across
    workers; results are collected as futures complete and then sorted by
    split_id before output, ensuring a deterministic CSV regardless of worker
    completion order.

Baseline comparison note
------------------------
    The seed Fisher baseline computed here (seeds as the query set, n = number
    of seeds in gene universe) is a different statistical object from the
    Analysis 4 / Table 4 Fisher baseline (graph neighborhood as query).
    Results are internally consistent across splits and suitable for reporting
    "BIFO vs Fisher within each split", but the absolute Fisher numbers must
    not be directly compared to Table 4 values without explanation.

Design notes
------------
    - GSEA baseline intentionally omitted from the resampling loop.
      The GSEA AP gap (BIFO ~0.40 vs raw GSEA ~0.12 on the primary split)
      is too large to be reversed by split variation; including GSEA would
      add ~2.3 s/split and inflate runtime ~30×.
    - All splits evaluated on the same frozen 550-pathway universe
      (via --bifo-scores), matching the primary benchmark exactly.
    - Supports claim: "The CHD discovery benchmark is not dependent on a
      single hand-selected 10/5 split of the 15-gene CHD pool."
"""

import argparse
import csv
import json
import math
import time
import numpy as np

# math.comb added in Python 3.8 (shim below); str.removesuffix replaced in clean()
if not hasattr(math, 'comb'):
    import functools, operator
    def _comb(n, k):
        if k < 0 or k > n:
            return 0
        k = min(k, n - k)
        num = functools.reduce(operator.mul, range(n, n - k, -1), 1)
        den = functools.reduce(operator.mul, range(1, k + 1), 1)
        return num // den
    math.comb = _comb


from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path
from scipy import sparse
from scipy.stats import hypergeom
from typing import Dict, List, Set, Tuple

# ─── CHD gene pool (15 genes) ─────────────────────────────────────────────────

CHD_POOL = [
    ('C1414995', 'GATA4'),
    ('C1413784', 'NKX2-5'),
    ('C1420615', 'TBX5'),
    ('C1334889', 'NOTCH1'),
    ('C1417767', 'NOTCH2'),
    ('C1415466', 'HAND1'),
    ('C1415467', 'HAND2'),
    ('C1417541', 'MYH6'),
    ('C1414996', 'GATA6'),
    ('C1420603', 'TBX1'),
    ('C1424490', 'ZFPM2'),
    ('C1417542', 'MYH7'),
    ('C1335280', 'PTPN11'),
    ('C1416525', 'JAG1'),
    ('C1333569', 'FLT4'),
]

# The original (primary benchmark) split — seeds are first 10, heldout last 5
ORIGINAL_SEEDS  = frozenset(c for c, _ in CHD_POOL[:10])
ORIGINAL_HELDOUT = frozenset(c for c, _ in CHD_POOL[10:])

POOL_CUIS = [c for c, _ in CHD_POOL]
POOL_SYMS = {c: g for c, g in CHD_POOL}

# ─── Utilities ────────────────────────────────────────────────────────────────

def clean(raw: str) -> str:
    s = raw.strip()
    return s[:-4] if s.endswith(' CUI') else s


def read_node_index(path: str) -> Dict[str, int]:
    return {clean(k): int(v) for k, v in json.load(open(path)).items()}


def row_normalize(A):
    """Row-normalize a sparse matrix; rows with zero sum left as zero."""
    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    D_inv = sparse.diags(1.0 / row_sums)
    return D_inv @ A


def ppr(A_tilde_T, seed_idx: List[int], n: int, alpha: float = 0.5,
        tol: float = 1e-10, max_iter: int = 500) -> np.ndarray:
    """Personalized PageRank. A_tilde_T = transpose of row-normalized A."""
    s = np.zeros(n, dtype=np.float64)
    if seed_idx:
        s[seed_idx] = 1.0 / len(seed_idx)
    f = s.copy()
    for _ in range(max_iter):
        nxt = (1.0 - alpha) * (A_tilde_T @ f) + alpha * s
        if np.linalg.norm(nxt - f, ord=1) < tol:
            return nxt
        f = nxt
    return f


def bh_correct(pvalues: List[float]) -> List[float]:
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed)):
        rank_from_end = n - rank
        adj = min(pval * n / rank_from_end, prev)
        prev = adj
        adjusted[orig_idx] = adj
    return adjusted


# ─── Build PPR operator ───────────────────────────────────────────────────────

def build_operator(kept_edges_path: str, node_index: Dict[str, int],
                   propagating_only: bool = True) -> sparse.csr_matrix:
    """Build sparse adjacency matrix from kept_edges.csv or edges_merged.csv.

    If propagating_only=True (default), only rows where the 'propagating'
    column is True are included. If the column is absent (e.g. edges_merged.csv),
    all rows are included.
    """
    n = len(node_index)
    rows_list, cols_list = [], []
    with open(kept_edges_path, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        has_prop = 'propagating' in (reader.fieldnames or [])
        for row in reader:
            if has_prop and propagating_only:
                prop = row.get('propagating', '').strip().lower()
                if prop not in ('true', '1', 'yes'):
                    continue
            src = clean(row.get('source', row.get('subject_cui', '')))
            tgt = clean(row.get('target', row.get('object_cui', '')))
            if src in node_index and tgt in node_index:
                rows_list.append(node_index[src])
                cols_list.append(node_index[tgt])
    data = np.ones(len(rows_list), dtype=np.float64)
    A = sparse.csr_matrix((data, (rows_list, cols_list)), shape=(n, n))
    return A


# ─── Build membership map ────────────────────────────────────────────────────

MEMBERSHIP_PREDS_FWD = {
    'pathway_associated_with_gene', 'has_signature_gene', 'process_involves_gene'
}
MEMBERSHIP_PREDS_REV = {
    'inverse_pathway_associated_with_gene', 'inverse_has_signature_gene',
    'gene_plays_role_in_process'
}


def build_membership(edges_merged_path: str, min_members: int, max_members: int,
                     excluded_patterns: List[str] = None) -> Tuple[
                         Dict[str, Set[str]], Dict[str, str], Dict[str, str], Set[str]]:
    """
    Returns:
        membership      : {pathway_id: frozenset of gene_ids}
        pathway_names   : {pathway_id: name}
        pathway_sabs    : {pathway_id: sab}
        gene_universe   : set of all C-prefixed gene nodes
    """
    raw_membership: Dict[str, Set[str]] = defaultdict(set)
    gene_universe: Set[str] = set()

    with open(edges_merged_path, encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            src = clean(row['source'])
            tgt = clean(row['target'])
            pred = row['predicate'].strip()
            if src.startswith('C'): gene_universe.add(src)
            if tgt.startswith('C'): gene_universe.add(tgt)
            if pred in MEMBERSHIP_PREDS_FWD:
                raw_membership[src].add(tgt)
            elif pred in MEMBERSHIP_PREDS_REV:
                raw_membership[tgt].add(src)

    # Filter by size and name patterns
    excluded = excluded_patterns or []
    membership = {}
    for pw, genes in raw_membership.items():
        if not (min_members <= len(genes) <= max_members):
            continue
        name = pathway_names.get(pw, pw) if 'pathway_names' in dir() else pw
        if any(pat in pw for pat in excluded):
            continue
        membership[pw] = frozenset(genes)

    return membership, gene_universe


def filter_to_bifo_universe(bifo_scores_path: str, membership: Dict) -> Dict:
    """Restrict membership to the BIFO-scored pathway universe."""
    bifo_ids = set()
    with open(bifo_scores_path) as f:
        for row in csv.DictReader(f):
            bifo_ids.add(row['concept_id'])
    return {k: v for k, v in membership.items() if k in bifo_ids}


# ─── Metrics ──────────────────────────────────────────────────────────────────

def ndcg_at_k(relevance: List[int], k: int) -> float:
    k = min(k, len(relevance))
    dcg = sum(r / math.log2(i + 2) for i, r in enumerate(relevance[:k]))
    n_rel = sum(relevance)
    ideal = sum(1.0 / math.log2(i + 2) for i in range(min(n_rel, k)))
    return dcg / ideal if ideal > 0 else 0.0


def average_precision(relevance: List[int]) -> float:
    hits, ap = 0, 0.0
    for i, r in enumerate(relevance):
        if r:
            hits += 1
            ap += hits / (i + 1)
    n_rel = sum(relevance)
    return ap / n_rel if n_rel > 0 else 0.0


def ranking_metrics(ranked_ids: List[str], ref: Set[str],
                    ks: Tuple[int, ...] = (10, 20, 50)) -> Dict:
    n = len(ranked_ids)
    n_ref = sum(1 for x in ranked_ids if x in ref)
    bg = n_ref / n if n > 0 else 0.0
    rel = [1 if x in ref else 0 for x in ranked_ids]
    m = {'n_scored': n, 'n_reference': n_ref, 'background_rate': bg}
    for k in ks:
        hits = sum(rel[:k])
        prec = hits / k
        m[f'p_at_{k}'] = prec
        m[f'enrich_at_{k}'] = prec / bg if bg > 0 else float('nan')
        m[f'recall_at_{k}'] = hits / n_ref if n_ref > 0 else 0.0
        m[f'ndcg_at_{k}'] = ndcg_at_k(rel, k)
    ref_ranks = [i + 1 for i, r in enumerate(rel) if r]
    m['mean_ref_rank'] = float(np.mean(ref_ranks)) if ref_ranks else float('nan')
    m['average_precision'] = average_precision(rel)
    return m


def score_pathways(f: np.ndarray, node_index: Dict[str, int],
                   membership: Dict[str, frozenset]) -> Dict[str, float]:
    """Compute degree_norm score: direct_ppr / sqrt(|members|)."""
    scores = {}
    for pw_id, members in membership.items():
        if pw_id not in node_index:
            scores[pw_id] = 0.0
            continue
        direct = float(f[node_index[pw_id]])
        n_members = len(members)
        scores[pw_id] = direct / math.sqrt(n_members) if n_members > 0 else 0.0
    return scores


def fisher_enrichment(query: Set[str], membership: Dict[str, frozenset],
                      gene_universe: Set[str]) -> Dict[str, float]:
    """Hypergeometric enrichment, returns {pw_id: pvalue}."""
    N = len(gene_universe)
    n = len(query & gene_universe)
    results = {}
    for pw_id, members in membership.items():
        K = len(members & gene_universe)
        k = len(members & query)
        results[pw_id] = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
    return results


def gsea_enrichment(gene_scores: Dict[str, float],
                    membership: Dict[str, frozenset],
                    gene_universe: Set[str]) -> Dict[str, float]:
    """Preranked GSEA enrichment score per pathway."""
    scored = [(g, s) for g, s in gene_scores.items() if g in gene_universe]
    scored.sort(key=lambda x: -x[1])
    ranked_genes = [g for g, _ in scored]
    ranked_scores = np.array([s for _, s in scored])
    N = len(ranked_genes)
    results = {}
    for pw_id, members in membership.items():
        pw_in_u = members & gene_universe
        if not pw_in_u:
            results[pw_id] = 0.0
            continue
        hits = np.array([1 if g in pw_in_u else 0 for g in ranked_genes])
        n_hit = hits.sum()
        if n_hit == 0:
            results[pw_id] = 0.0
            continue
        hw = np.where(hits == 1, np.abs(ranked_scores), 0.0)
        hw_sum = hw.sum()
        if hw_sum == 0:
            hw = hits.astype(float); hw_sum = float(n_hit)
        miss_pen = 1.0 / (N - n_hit) if N > n_hit else 0.0
        rs, max_dev = 0.0, 0.0
        for i in range(N):
            rs += hw[i] / hw_sum if hits[i] else -miss_pen
            if abs(rs) > abs(max_dev): max_dev = rs
        results[pw_id] = float(max_dev)
    return results


def auprc(f: np.ndarray, heldout_idx: List[int],
          seed_idx: List[int]) -> float:
    """AUPRC for held-out gene recovery."""
    if not heldout_idx:
        return float('nan')
    exclude = set(seed_idx) | set(heldout_idx)
    labels, scores_out = [], []
    for i in range(len(f)):
        if i in exclude:
            continue
        labels.append(0)
        scores_out.append(float(f[i]))
    for i in heldout_idx:
        labels.append(1)
        scores_out.append(float(f[i]))
    # Simple AUPRC via sorting
    paired = sorted(zip(scores_out, labels), reverse=True)
    tp, fp, prev_recall = 0, 0, 0.0
    n_pos = sum(labels)
    ap = 0.0
    for _, label in paired:
        if label: tp += 1
        else: fp += 1
        recall = tp / n_pos if n_pos > 0 else 0.0
        prec = tp / (tp + fp)
        ap += prec * (recall - prev_recall)
        prev_recall = recall
    return ap


# ─── Parallel worker state ─────────────────────────────────────────────────────
# These globals are populated by _worker_init() in each spawned process.

_W_A_T      = None  # conditioned operator transpose (CSR)
_W_A_RAW_T  = None  # raw operator transpose (CSR)
_W_N        = None  # number of nodes
_W_MEM      = None  # membership dict
_W_CHD      = None  # CHD reference frozenset
_W_GU       = None  # gene universe frozenset
_W_NI       = None  # node_index dict
_W_POOL     = None  # CHD pool CUI list (15 genes)
_W_ORIG     = None  # original split seed frozenset (indices)
_W_ALPHA    = None  # PPR restart probability


def _worker_init(at_data, at_indices, at_indptr, at_shape,
                 araw_data, araw_indices, araw_indptr, araw_shape,
                 n_nodes, membership_items, chd_ref_list,
                 gene_universe_list, node_index_items,
                 pool_cuis, orig_seed_indices, alpha):
    """Initializer for each worker process. Reconstructs all shared state."""
    global _W_A_T, _W_A_RAW_T, _W_N, _W_MEM, _W_CHD, _W_GU
    global _W_NI, _W_POOL, _W_ORIG, _W_ALPHA
    _W_A_T     = sparse.csr_matrix((at_data,   at_indices,   at_indptr),   shape=at_shape)
    _W_A_RAW_T = sparse.csr_matrix((araw_data, araw_indices, araw_indptr), shape=araw_shape)
    _W_N       = n_nodes
    _W_MEM     = {k: frozenset(v) for k, v in membership_items}
    _W_CHD     = frozenset(chd_ref_list)
    _W_GU      = frozenset(gene_universe_list)
    _W_NI      = dict(node_index_items)
    _W_POOL    = pool_cuis
    _W_ORIG    = frozenset(orig_seed_indices)
    _W_ALPHA   = alpha


def _process_batch(combo_batch: List[Tuple[int, ...]]) -> List[dict]:
    """
    Process a batch of split combos in a worker process.
    combo_batch: list of tuples of pool indices (e.g. (0,1,2,...,9) for original split).
    Returns a list of result dicts (one per split).
    """
    results = []
    for combo in combo_batch:
        seed_set   = frozenset(combo)
        heldout_set = frozenset(range(15)) - seed_set
        seed_cuis   = [_W_POOL[i] for i in sorted(seed_set)]
        heldout_cuis = [_W_POOL[i] for i in sorted(heldout_set)]

        s_idx = [_W_NI[c] for c in seed_cuis   if c in _W_NI]
        h_idx = [_W_NI[c] for c in heldout_cuis if c in _W_NI]
        if not s_idx:
            continue

        f_c = ppr(_W_A_T,     s_idx, _W_N, alpha=_W_ALPHA)
        f_r = ppr(_W_A_RAW_T, s_idx, _W_N, alpha=_W_ALPHA)

        # ── Pathway scoring ──
        pw_c = score_pathways(f_c, _W_NI, _W_MEM)
        pw_r = {pw: float(f_r[_W_NI[pw]]) if pw in _W_NI else 0.0 for pw in _W_MEM}

        ranked_c = sorted(pw_c, key=pw_c.get, reverse=True)
        ranked_r = sorted(pw_r, key=pw_r.get, reverse=True)
        mc = ranking_metrics(ranked_c, _W_CHD)
        mr = ranking_metrics(ranked_r, _W_CHD)
        rank_imp = mr['mean_ref_rank'] - mc['mean_ref_rank']

        # ── Seed Fisher baseline ──
        sc = frozenset(seed_cuis)
        N  = len(_W_GU)
        n_q = len(sc & _W_GU)
        sf_pvals = {}
        for pw_id, members in _W_MEM.items():
            K = len(members & _W_GU)
            k = len(members & sc)
            sf_pvals[pw_id] = hypergeom.sf(k - 1, N, K, n_q) if k > 0 else 1.0
        sf_ranked = sorted(sf_pvals, key=sf_pvals.get)
        sf_m = ranking_metrics(sf_ranked, _W_CHD)

        # ── Gene-level AUPRC ──
        gene_ap = auprc(f_c, h_idx, s_idx)

        # split_id: lexicographically stable identifier built from sorted pool
        # indices (immune to CUI formatting). Used to sort output CSV deterministically.
        split_id = '-'.join(str(i) for i in sorted(seed_set))

        results.append({
            'split_id':         split_id,
            'is_original':      int(seed_set == _W_ORIG),
            'seeds':            ','.join(_W_POOL[i] for i in sorted(seed_set)),
            'heldout':          ','.join(_W_POOL[i] for i in sorted(heldout_set)),
            'bifo_p10':         mc['p_at_10'],
            'bifo_p20':         mc['p_at_20'],
            'bifo_enrich10':    mc['enrich_at_10'],
            'bifo_ndcg10':      mc['ndcg_at_10'],
            'bifo_ap':          mc['average_precision'],
            'bifo_mean_rank':   mc['mean_ref_rank'],
            'rank_improvement': rank_imp,
            'gene_auprc':       gene_ap,
            'sf_p10':           sf_m['p_at_10'],
            'sf_ap':            sf_m['average_precision'],
            'bifo_beats_sf_ap': int(mc['average_precision'] > sf_m['average_precision']),
        })
    return results


# ─── Main ─────────────────────────────────────────────────────────────────────


# ─── Main ─────────────────────────────────────────────────────────────────────

def run(args):
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed

    print("=" * 70)
    print("CHD Exhaustive Resampling Analysis")
    print(f"  Pool: {len(CHD_POOL)} genes  |  "
          f"C({len(CHD_POOL)},10) = {math.comb(len(CHD_POOL),10):,} splits")
    print("=" * 70)

    # ── Build operators ──────────────────────────────────────────────────────
    print("\n[1/5] Loading node index and building conditioned PPR operator...")
    node_index = read_node_index(args.node_index)
    n = len(node_index)
    A      = build_operator(args.kept_edges, node_index, propagating_only=True)
    A_T    = row_normalize(A).T.tocsr()
    print(f"      {A.nnz:,} propagating edges, {n:,} nodes")

    print("[2/5] Building raw (full) operator from edges_merged...")
    A_raw  = build_operator(args.edges_merged, node_index, propagating_only=False)
    A_raw_T = row_normalize(A_raw).T.tocsr()
    print(f"      {A_raw.nnz:,} edges")

    # ── Pathway data ─────────────────────────────────────────────────────────
    print("[3/5] Building pathway membership map...")
    _excluded_patterns = ['_Q2', '_Q3', '_Q4', '_Q5', '_Q6', 'MIR']
    membership, gene_universe = build_membership(
        args.edges_merged, args.min_members, args.max_members,
        excluded_patterns=_excluded_patterns)
    if args.bifo_scores and Path(args.bifo_scores).exists():
        membership = filter_to_bifo_universe(args.bifo_scores, membership)
    print(f"      {len(membership):,} pathways (frozen BIFO universe)")

    chd_ref: Set[str] = set()
    for line in Path(args.chd_pathways).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith('#'):
            chd_ref.add(clean(line.split()[0]))
    chd_ref &= set(membership.keys())
    print(f"      CHD reference: {len(chd_ref)} pathways "
          f"(background {len(chd_ref)/len(membership):.3f})")

    if any(c not in node_index for c, _ in CHD_POOL):
        missing = [c for c, _ in CHD_POOL if c not in node_index]
        print(f"  WARNING: {len(missing)} pool genes not in node index: {missing}")

    # Gene universe: nodes with edges, not all node_index entries.
    # Matches baseline_enrichment.py universe construction.
    gene_universe_set = gene_universe  # built from edge endpoints in build_membership

    # ── Enumerate splits ──────────────────────────────────────────────────────
    all_cuis   = [c for c, _ in CHD_POOL]
    all_splits = list(combinations(range(len(all_cuis)), 10))
    if args.max_splits:
        all_splits = all_splits[:args.max_splits]

    n_splits  = len(all_splits)
    n_workers = args.n_cores or os.cpu_count() or 1

    print(f"\n[4/5] Running {n_splits:,} splits on {n_workers} worker(s)...")

    # Serialise operators as CSR component arrays (the only way to
    # pass scipy sparse matrices safely to subprocess workers).
    init_args = (
        A_T.data,     A_T.indices,     A_T.indptr,     A_T.shape,
        A_raw_T.data, A_raw_T.indices, A_raw_T.indptr, A_raw_T.shape,
        n,
        list(membership.items()),          # [(pw_id, frozenset), ...]
        list(chd_ref),                     # [cui, ...]
        list(gene_universe_set),           # [cui, ...]
        list(node_index.items()),          # [(cui, idx), ...]
        all_cuis,                          # [cui, ...]  len=15
        list(range(10)),                   # original split = first 10
        args.alpha,
    )

    # Chunk splits evenly across workers
    chunk_size = max(1, math.ceil(n_splits / n_workers))
    chunks = [all_splits[i:i+chunk_size]
              for i in range(0, n_splits, chunk_size)]

    results_rows: List[dict] = []
    t0 = time.time()

    def _collect(batch_result):
        results_rows.extend(batch_result)
        done = len(results_rows)
        elapsed = time.time() - t0
        eta = elapsed / done * (n_splits - done) if done < n_splits else 0
        print(f"      {done:4d}/{n_splits}  "
              f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")

    if n_workers == 1:
        # Serial path: avoids subprocess overhead, easier to debug
        _worker_init(*init_args)
        for chunk in chunks:
            _collect(_process_batch(chunk))
    else:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_worker_init,
            initargs=init_args,
        ) as ex:
            for future in as_completed(
                ex.submit(_process_batch, chunk) for chunk in chunks
            ):
                _collect(future.result())

    elapsed = time.time() - t0
    n_done  = len(results_rows)
    print(f"      Done. {n_done} splits in {elapsed:.1f}s "
          f"({elapsed/n_done*1000:.0f} ms/split avg)")

    # ── Aggregate ─────────────────────────────────────────────────────────────
    print("\n[5/5] Aggregating and writing results...")

    def stats(key: str) -> dict:
        vals = [r[key] for r in results_rows
                if isinstance(r.get(key), (int, float))
                and not math.isnan(r[key])]
        if not vals:
            return {}
        a = np.array(vals, dtype=float)
        return {
            'n':      len(a),
            'mean':   float(a.mean()),
            'sd':     float(a.std(ddof=1)) if len(a) > 1 else 0.0,
            'min':    float(a.min()),
            'p25':    float(np.percentile(a, 25)),
            'median': float(np.median(a)),
            'p75':    float(np.percentile(a, 75)),
            'max':    float(a.max()),
        }

    numeric_keys = [k for k in results_rows[0]
                    if k not in ('is_original', 'seeds', 'heldout')]

    # Sort by split_id for deterministic CSV output regardless of worker
    # completion order.  Must happen before orig_row extraction and CSV write.
    results_rows.sort(key=lambda r: r['split_id'])

    orig_row = next((r for r in results_rows if r['is_original']), None)

    summary = {
        'n_splits_total':  math.comb(len(CHD_POOL), 10),
        'n_splits_run':    n_done,
        'exhaustive':      n_done == math.comb(len(CHD_POOL), 10),
        'alpha':           args.alpha,
        'pathway_universe':len(membership),
        'chd_reference':   len(chd_ref),
        'background_rate': len(chd_ref) / len(membership),
        'original_split':  orig_row,
        'robustness': {
            'bifo_p10_ge_0.5':   sum(1 for r in results_rows if r['bifo_p10'] >= 0.5),
            'bifo_p10_ge_0.3':   sum(1 for r in results_rows if r['bifo_p10'] >= 0.3),
            'rank_imp_positive': sum(1 for r in results_rows if r['rank_improvement'] > 0),
            'bifo_beats_sf_ap':  sum(r['bifo_beats_sf_ap'] for r in results_rows),
        },
        'metrics': {k: stats(k) for k in numeric_keys},
    }

    # Write per-split CSV
    csv_path = Path(args.out_csv)
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results_rows[0].keys()))
        writer.writeheader()
        writer.writerows(results_rows)

    # Write summary JSON
    Path(args.out_json).write_text(json.dumps(summary, indent=2))

    # ── Print summary ─────────────────────────────────────────────────────────
    print()
    print("=" * 72)
    print(f"RESAMPLING SUMMARY  "
          f"({'exhaustive' if summary['exhaustive'] else str(n_done)+' sampled'}, "
          f"{n_done} splits)")
    print("=" * 72)

    hdr = f"  {'Metric':40s}  {'Mean':>6s}  {'SD':>6s}  {'P25':>5s}  {'Median':>6s}  {'P75':>5s}"
    print(hdr)
    print("  " + "-" * 70)
    for key, label in [
        ('bifo_p10',         'BIFO P@10'),
        ('bifo_enrich10',    'BIFO Enrichment@10'),
        ('bifo_ndcg10',      'BIFO NDCG@10'),
        ('bifo_ap',          'BIFO Avg. Precision'),
        ('rank_improvement', 'BIFO Rank Improvement'),
        ('bifo_mean_rank',   'BIFO Mean CHD Rank'),
        ('gene_auprc',       'Conditioned AUPRC (gene)'),
        ('sf_p10',           'Seed Fisher P@10'),
        ('sf_ap',            'Seed Fisher AP'),
    ]:
        s = summary['metrics'].get(key, {})
        if s:
            print(f"  {label:40s}  {s['mean']:>6.3f}  {s['sd']:>6.3f}  "
                  f"{s['p25']:>5.3f}  {s['median']:>6.3f}  {s['p75']:>5.3f}")

    if orig_row:
        print()
        print("  Original split (primary benchmark):")
        for key, label in [
            ('bifo_p10',         'BIFO P@10'),
            ('bifo_ap',          'BIFO Avg. Precision'),
            ('rank_improvement', 'BIFO Rank Improvement'),
            ('gene_auprc',       'Conditioned AUPRC (gene)'),
        ]:
            v = orig_row.get(key, float('nan'))
            vals = sorted(r[key] for r in results_rows
                          if not math.isnan(r.get(key, float('nan'))))
            pct = sum(1 for x in vals if x <= v) / len(vals) * 100 if vals else float('nan')
            print(f"    {label:40s}  {v:>6.3f}  (pct rank: {pct:.0f}%)")

    rob = summary['robustness']
    n = n_done
    print()
    print(f"  Robustness ({n} splits):")
    print(f"    BIFO P@10 >= 0.50:           {rob['bifo_p10_ge_0.5']:4d}/{n} "
          f"({100*rob['bifo_p10_ge_0.5']/n:.1f}%)")
    print(f"    BIFO P@10 >= 0.30:           {rob['bifo_p10_ge_0.3']:4d}/{n} "
          f"({100*rob['bifo_p10_ge_0.3']/n:.1f}%)")
    print(f"    Rank improvement > 0:        {rob['rank_imp_positive']:4d}/{n} "
          f"({100*rob['rank_imp_positive']/n:.1f}%)")
    print(f"    BIFO beats Fisher (AP):      {rob['bifo_beats_sf_ap']:4d}/{n} "
          f"({100*rob['bifo_beats_sf_ap']/n:.1f}%)")
    print(f"    NOTE: Resampling Fisher ≠ Analysis 4 Table 4 Fisher.")
    print(f"    Resampling uses seeds as query; Analysis 4 uses graph neighborhood.")
    print()
    print(f"  Output: {csv_path} ({csv_path.stat().st_size//1024} KB)")
    print(f"          {args.out_json}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="CHD exhaustive resampling — all C(15,10)=3,003 splits, "
                    "in-memory, optional parallel execution via --n-cores"
    )
    p.add_argument('--kept-edges',   required=True,
                   help='results_kept_edges.csv from bifo_conditioning.py')
    p.add_argument('--edges-merged', required=True,
                   help='edges_merged.csv (used for raw arm + membership)')
    p.add_argument('--node-index',   required=True,
                   help='results_node_index.json from bifo_conditioning.py')
    p.add_argument('--bifo-scores',  default=None,
                   help='pathway_scores_full.csv — restricts to frozen BIFO universe')
    p.add_argument('--chd-pathways', required=True,
                   help='CHD reference pathway CUI list')
    p.add_argument('--out-csv',      required=True,
                   help='Per-split metrics CSV output')
    p.add_argument('--out-json',     required=True,
                   help='Aggregated summary JSON output')
    p.add_argument('--alpha',        type=float, default=0.5)
    p.add_argument('--min-members',  type=int,   default=8)
    p.add_argument('--max-members',  type=int,   default=300)
    p.add_argument('--max-splits',   type=int,   default=None,
                   help='Limit to first N splits for testing (default: all 3,003)')
    p.add_argument('--n-cores',      type=int,   default=None,
                   help='Parallel workers (default: all CPU cores). '
                        'Use 1 for serial/debug mode.')
    run(p.parse_args())


if __name__ == '__main__':
    main()
