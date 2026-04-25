#!/usr/bin/env python3
"""
kf_resampling.py

Bootstrap resampling for Kids First cohort BIFO analyses.

Unlike the CHD curated benchmark (exhaustive C(15,10)=3,003 splits),
the KF cohort pools are too large for exhaustive enumeration. Instead
we draw B random subsamples at each of several seed sizes, run the
full BIFO pipeline on each, and report the distribution of pathway
metrics across draws.

This answers: "Is pathway prioritization stable across different subsets
of the variant-derived gene list, and across different seed set sizes?"

Usage
-----
  python kf_resampling.py \\
    --kept-edges   kf_chd_results/results_kept_edges.csv \\
    --edges-merged kf_chd_edges_all.csv \\
    --node-index   kf_chd_results/results_node_index.json \\
    --bifo-scores  kf_chd_results/pathway_scores_standard.csv \\
    --seed-pool    kf_chd_seeds.txt \\
    --ref-pathways kf_chd_cilia_reference.txt \\
    --out-csv      kf_chd_results/resampling_results.csv \\
    --out-json     kf_chd_results/resampling_summary.json \\
    --seed-sizes   10 20 30 \\
    --n-boots      500 \\
    --n-cores      192 \\
    --primary-size 20
"""

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s


# Python 3.7 compatibility: str.removesuffix was added in 3.9
def _removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

import csv
import json
import math
import os
import random
import time
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from scipy import sparse
from scipy.stats import hypergeom
from typing import Dict, List, Set, Tuple

# ─── Compatibility shims ──────────────────────────────────────────────────────
if not hasattr(math, 'comb'):
    import functools, operator
    def _comb(n, k):
        if k < 0 or k > n: return 0
        k = min(k, n - k)
        return (functools.reduce(operator.mul, range(n, n-k, -1), 1) //
                functools.reduce(operator.mul, range(1, k+1), 1))
    math.comb = _comb

if not hasattr(str, 'removesuffix'):
    def _removesuffix(s, suffix):
        return s[:-len(suffix)] if suffix and s.endswith(suffix) else s


# ─── Graph helpers ────────────────────────────────────────────────────────────

def clean(raw: str) -> str:
    return _removesuffix(raw.strip(), ' CUI')


def read_node_index(path: str) -> Dict[str, int]:
    return {clean(k): int(v) for k, v in json.load(open(path)).items()}


def build_operator(edges_csv: str, node_index: Dict[str, int],
                   propagating_only: bool = True) -> sparse.csr_matrix:
    n = len(node_index)
    rows, cols = [], []
    with open(edges_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if propagating_only:
                prop = row.get('propagating', '').strip().lower()
                if prop not in ('true', '1', 'yes'):
                    continue
            s = clean(row.get('source', row.get('subject', '')))
            t = clean(row.get('target', row.get('object', '')))
            if s in node_index and t in node_index:
                rows.append(node_index[s])
                cols.append(node_index[t])
    data = np.ones(len(rows), dtype=np.float32)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n))


def row_normalize(A: sparse.csr_matrix) -> sparse.csr_matrix:
    row_sums = np.array(A.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1.0
    inv = sparse.diags(1.0 / row_sums)
    return inv @ A


def ppr(A_T: sparse.csr_matrix, seed_idx: List[int], n: int,
        alpha: float = 0.5, tol: float = 1e-10, max_iter: int = 500) -> np.ndarray:
    s = np.zeros(n, dtype=np.float64)
    if not seed_idx:
        return s
    for i in seed_idx:
        s[i] = 1.0 / len(seed_idx)
    f = s.copy()
    for _ in range(max_iter):
        f_new = (1 - alpha) * A_T.dot(f) + alpha * s
        if np.abs(f_new - f).sum() < tol:
            return f_new
        f = f_new
    return f


def build_membership(edges_csv: str, min_members: int = 8,
                     max_members: int = 300,
                     excluded_patterns: list = None):
    """
    Build pathway membership map from edges file.

    Applies size filters (min_members, max_members) and optionally excludes
    pathways whose IDs match any string in excluded_patterns. This is used to
    remove gene expression quantile sets (e.g. _Q2 through _Q6) and microRNA
    target sets (MIR) from the resampling universe, as these represent
    statistical partitions rather than curated biological programs. See
    Methods §5.4 and §9 for the name-pattern filter rationale.
    """
    MEMBERSHIP_PREDS = {
        'pathway_associated_with_gene', 'inverse_pathway_associated_with_gene',
        'has_signature_gene', 'inverse_has_signature_gene',
        'process_involves_gene', 'gene_plays_role_in_process',
    }
    pw_to_genes: Dict[str, Set[str]] = defaultdict(set)
    gene_universe: Set[str] = set()
    with open(edges_csv) as f:
        for row in csv.DictReader(f):
            pred = row.get('predicate', '').strip()
            if pred not in MEMBERSHIP_PREDS:
                continue
            src = clean(row.get('source', row.get('subject', '')))
            tgt = clean(row.get('target', row.get('object', '')))
            if pred in ('pathway_associated_with_gene', 'has_signature_gene',
                        'process_involves_gene'):
                pw_to_genes[src].add(tgt)
                gene_universe.add(tgt)
            else:
                pw_to_genes[tgt].add(src)
                gene_universe.add(src)
    excl = excluded_patterns or []
    membership = {pw: frozenset(genes)
                  for pw, genes in pw_to_genes.items()
                  if min_members <= len(genes) <= max_members
                  and not any(pat in pw for pat in excl)}
    return membership, gene_universe


def filter_to_bifo_universe(bifo_scores_csv: str,
                             membership: Dict) -> Dict:
    scored = set()
    with open(bifo_scores_csv) as f:
        for row in csv.DictReader(f):
            cui = (row.get('concept_id') or row.get('pathway_cui') or
                   row.get('CUI') or row.get('cui') or '')
            if cui:
                scored.add(clean(cui))
    return {pw: members for pw, members in membership.items() if pw in scored}


def score_pathways(f: np.ndarray, node_index: Dict[str, int],
                   membership: Dict) -> Dict[str, float]:
    scores = {}
    for pw, members in membership.items():
        if pw in node_index:
            n_members = len(members)
            scores[pw] = float(f[node_index[pw]]) / math.sqrt(max(n_members, 1))
        else:
            scores[pw] = 0.0
    return scores


def ranking_metrics(ranked: List[str], ref: Set[str],
                    bg_rate: float = None) -> Dict:
    n_ref = len(ref)
    n_uni = len(ranked)
    bg = bg_rate or (n_ref / n_uni if n_uni > 0 else 0)

    hits_10 = sum(1 for p in ranked[:10] if p in ref)
    p10 = hits_10 / 10
    enr10 = p10 / bg if bg > 0 else 0

    # NDCG@10
    dcg = sum((1 if ranked[i] in ref else 0) / math.log2(i + 2)
               for i in range(min(10, len(ranked))))
    ideal = sum(1 / math.log2(i + 2) for i in range(min(n_ref, 10)))
    ndcg10 = dcg / ideal if ideal > 0 else 0

    # AP and mean rank
    hits, ap_sum = 0, 0.0
    ref_ranks = []
    for i, p in enumerate(ranked):
        if p in ref:
            hits += 1
            ap_sum += hits / (i + 1)
            ref_ranks.append(i + 1)
    ap = ap_sum / n_ref if n_ref > 0 else 0
    mean_rank = sum(ref_ranks) / len(ref_ranks) if ref_ranks else n_uni + 1

    return {'p_at_10': p10, 'enrich_at_10': enr10, 'ndcg_at_10': ndcg10,
            'average_precision': ap, 'mean_ref_rank': mean_rank}


# ─── Worker globals ───────────────────────────────────────────────────────────

_W: dict = {}


def _worker_init(at_data, at_idx, at_ptr, at_shape,
                 araw_data, araw_idx, araw_ptr, araw_shape,
                 n, membership_items, ref_list, gene_universe_list,
                 node_index_items, alpha):
    global _W
    _W['A_T']     = sparse.csr_matrix((at_data,   at_idx,   at_ptr),   shape=at_shape)
    _W['A_raw_T'] = sparse.csr_matrix((araw_data, araw_idx, araw_ptr), shape=araw_shape)
    _W['n']       = n
    _W['mem']     = {k: frozenset(v) for k, v in membership_items}
    _W['ref']     = frozenset(ref_list)
    _W['gu']      = frozenset(gene_universe_list)
    _W['ni']      = dict(node_index_items)
    _W['alpha']   = alpha


def _process_batch(tasks: List[Tuple[int, List[str], int]]) -> List[dict]:
    """
    tasks: list of (boot_id, seed_cui_list, seed_size)
    """
    results = []
    for boot_id, seed_cuis, seed_size in tasks:
        s_idx = [_W['ni'][c] for c in seed_cuis if c in _W['ni']]
        if not s_idx:
            continue

        f_c = ppr(_W['A_T'],     s_idx, _W['n'], _W['alpha'])
        f_r = ppr(_W['A_raw_T'], s_idx, _W['n'], _W['alpha'])

        pw_c = score_pathways(f_c, _W['ni'], _W['mem'])
        pw_r = {pw: float(f_r[_W['ni'][pw]]) if pw in _W['ni'] else 0.0
                for pw in _W['mem']}

        ranked_c = sorted(pw_c, key=pw_c.get, reverse=True)
        ranked_r = sorted(pw_r, key=pw_r.get, reverse=True)
        mc = ranking_metrics(ranked_c, _W['ref'])
        mr = ranking_metrics(ranked_r, _W['ref'])
        rank_imp = mr['mean_ref_rank'] - mc['mean_ref_rank']

        # Seed Fisher
        sc = frozenset(seed_cuis)
        N = len(_W['gu'])
        n_q = len(sc & _W['gu'])
        sf_pvals = {}
        for pw_id, members in _W['mem'].items():
            K = len(members & _W['gu'])
            k = len(members & sc)
            sf_pvals[pw_id] = hypergeom.sf(k-1, N, K, n_q) if k > 0 else 1.0
        sf_ranked = sorted(sf_pvals, key=sf_pvals.get)
        sf_m = ranking_metrics(sf_ranked, _W['ref'])

        results.append({
            'boot_id':          boot_id,
            'seed_size':        seed_size,
            'seeds':            ','.join(seed_cuis),
            'bifo_p10':         mc['p_at_10'],
            'bifo_enrich10':    mc['enrich_at_10'],
            'bifo_ndcg10':      mc['ndcg_at_10'],
            'bifo_ap':          mc['average_precision'],
            'bifo_mean_rank':   mc['mean_ref_rank'],
            'rank_improvement': rank_imp,
            'sf_p10':           sf_m['p_at_10'],
            'sf_ap':            sf_m['average_precision'],
            'bifo_beats_sf_ap': int(mc['average_precision'] > sf_m['average_precision']),
        })
    return results


# ─── Main ─────────────────────────────────────────────────────────────────────

def run(args):
    print("=" * 70)
    print(f"KF Cohort Bootstrap Resampling")
    print(f"  Seed sizes: {args.seed_sizes}  |  Bootstraps: {args.n_boots}")
    print("=" * 70)

    # Load seed pool
    pool = []  # [(symbol, cui, n_carriers)]
    with open(args.seed_pool) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            pool.append(parts[0].strip())
    print(f"\nSeed pool: {len(pool)} genes")

    # Map symbols to CUIs via node index
    print("\n[1/5] Loading graph structures...")
    node_index = read_node_index(args.node_index)
    n = len(node_index)

    # Map seed pool symbols to CUIs
    # node_index keys are CUIs; we need to find CUIs for pool symbols
    # Try loading from a separate seed_cuis file if present
    pool_cuis = []
    cui_file = Path(args.seed_pool).with_suffix('').parent / \
               (Path(args.seed_pool).stem.replace('seeds', 'seed_cuis') + '.txt')
    if cui_file.exists():
        with open(cui_file) as f:
            for line in f:
                if line.startswith('#') or not line.strip(): continue
                cui = line.strip().split()[0].strip()
                pool_cuis.append(cui)
        print(f"  Loaded {len(pool_cuis)} seed CUIs from {cui_file.name}")
    else:
        # Fall back: filter node_index to seed symbols via PREF_TERM name pattern
        # node names in index are CUIs, so we need to cross-reference
        print(f"  WARNING: no seed CUI file found at {cui_file}")
        print(f"  Pass --seed-cuis explicitly or ensure seed_cuis.txt exists")
        return

    if len(pool_cuis) < min(args.seed_sizes):
        print(f"ERROR: only {len(pool_cuis)} CUIs but smallest seed size is {min(args.seed_sizes)}")
        return

    # Build operators
    A     = build_operator(args.kept_edges, node_index, propagating_only=True)
    A_T   = row_normalize(A).T.tocsr()
    print(f"  Conditioned: {A.nnz:,} propagating edges, {n:,} nodes")

    A_raw  = build_operator(args.edges_merged, node_index, propagating_only=False)
    A_raw_T = row_normalize(A_raw).T.tocsr()
    print(f"  Raw: {A_raw.nnz:,} edges")

    print("\n[2/5] Building pathway membership...")
    # Exclude quantile-binned gene sets and miRNA target sets from the
    # resampling universe — these are statistical partitions, not biological
    # programs. Consistent with name-pattern filter in score_pathways.py
    # and documented in Methods §5.4 and §9.
    _excluded_patterns = ['_Q2', '_Q3', '_Q4', '_Q5', '_Q6', 'MIR']
    membership, gene_universe = build_membership(args.edges_merged,
                                                  args.min_members, args.max_members,
                                                  excluded_patterns=_excluded_patterns)
    if args.bifo_scores and Path(args.bifo_scores).exists():
        membership = filter_to_bifo_universe(args.bifo_scores, membership)
    print(f"  {len(membership):,} pathways in universe")

    # Load reference set
    ref: Set[str] = set()
    with open(args.ref_pathways) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ref.add(clean(line.split()[0]))
    ref &= set(membership.keys())
    print(f"  Reference: {len(ref)} pathways  (background {len(ref)/len(membership):.3f})")

    # Build bootstrap tasks
    print(f"\n[3/5] Building {args.n_boots} bootstrap draws per seed size...")
    rng = random.Random(args.random_seed)
    tasks = []
    boot_id = 0
    for size in args.seed_sizes:
        if size > len(pool_cuis):
            print(f"  Skipping size={size} (pool only has {len(pool_cuis)} CUIs)")
            continue
        for _ in range(args.n_boots):
            sample = rng.sample(pool_cuis, size)
            tasks.append((boot_id, sample, size))
            boot_id += 1
    print(f"  Total tasks: {len(tasks)}")

    # Primary run (full seed set)
    primary_cuis = pool_cuis
    tasks.append((-1, primary_cuis, len(primary_cuis)))
    print(f"  + 1 primary run ({len(primary_cuis)} seeds)")

    # Worker init args
    n_workers = args.n_cores or os.cpu_count() or 1
    init_args = (
        A_T.data,     A_T.indices,     A_T.indptr,     A_T.shape,
        A_raw_T.data, A_raw_T.indices, A_raw_T.indptr, A_raw_T.shape,
        n,
        list(membership.items()),
        list(ref),
        list(gene_universe),
        list(node_index.items()),
        args.alpha,
    )

    chunk_size = max(1, math.ceil(len(tasks) / n_workers))
    chunks = [tasks[i:i+chunk_size] for i in range(0, len(tasks), chunk_size)]

    print(f"\n[4/5] Running {len(tasks):,} tasks on {n_workers} worker(s)...")
    results_rows = []
    t0 = time.time()

    if n_workers == 1:
        _worker_init(*init_args)
        for chunk in chunks:
            results_rows.extend(_process_batch(chunk))
    else:
        with ProcessPoolExecutor(max_workers=n_workers,
                                  initializer=_worker_init,
                                  initargs=init_args) as ex:
            for future in as_completed(ex.submit(_process_batch, c) for c in chunks):
                results_rows.extend(future.result())

    elapsed = time.time() - t0
    print(f"  Done. {len(results_rows)} runs in {elapsed:.1f}s")

    # Separate primary from bootstrap
    primary_row = next((r for r in results_rows if r['boot_id'] == -1), None)
    boot_rows   = [r for r in results_rows if r['boot_id'] != -1]

    print("\n[5/5] Aggregating results...")

    def stats(vals):
        if not vals: return {}
        a = np.array(vals, dtype=float)
        return {'mean': float(a.mean()), 'sd': float(a.std(ddof=1)) if len(a)>1 else 0.0,
                'min': float(a.min()), 'p25': float(np.percentile(a,25)),
                'median': float(np.median(a)), 'p75': float(np.percentile(a,75)),
                'max': float(a.max()), 'n': len(a)}

    by_size = defaultdict(list)
    for r in boot_rows:
        by_size[r['seed_size']].append(r)

    summary = {
        'cohort':          Path(args.seed_pool).stem.replace('_seeds',''),
        'n_boots':         args.n_boots,
        'seed_sizes':      args.seed_sizes,
        'pathway_universe':len(membership),
        'reference_size':  len(ref),
        'background_rate': len(ref)/len(membership) if membership else 0,
        'alpha':           args.alpha,
        'primary_run':     primary_row,
        'by_seed_size':    {},
        'overall': {
            'rank_imp_positive': sum(1 for r in boot_rows if r['rank_improvement']>0),
            'rank_imp_positive_pct': 100*sum(1 for r in boot_rows if r['rank_improvement']>0)/len(boot_rows) if boot_rows else 0,
            'bifo_beats_sf_pct': 100*sum(r['bifo_beats_sf_ap'] for r in boot_rows)/len(boot_rows) if boot_rows else 0,
        }
    }

    for size, rows in sorted(by_size.items()):
        summary['by_seed_size'][size] = {
            'n': len(rows),
            'bifo_p10':         stats([r['bifo_p10'] for r in rows]),
            'bifo_ap':          stats([r['bifo_ap'] for r in rows]),
            'rank_improvement': stats([r['rank_improvement'] for r in rows]),
            'sf_ap':            stats([r['sf_ap'] for r in rows]),
            'rank_imp_positive': sum(1 for r in rows if r['rank_improvement']>0),
            'bifo_beats_sf':     sum(r['bifo_beats_sf_ap'] for r in rows),
        }

    # Write outputs
    with open(args.out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(results_rows[0].keys()))
        writer.writeheader()
        writer.writerows(results_rows)

    Path(args.out_json).write_text(json.dumps(summary, indent=2))

    # Print summary
    print()
    print("=" * 70)
    print(f"BOOTSTRAP SUMMARY — {summary['cohort']}")
    print("=" * 70)
    for size, s in summary['by_seed_size'].items():
        ri = s['rank_improvement']
        ap = s['bifo_ap']
        print(f"\n  Seed size = {size}  (n={s['n']} bootstraps)")
        print(f"  {'Metric':25s}  {'Mean':>6s}  {'SD':>6s}  {'P25':>5s}  {'Median':>6s}  {'P75':>5s}")
        print(f"  {'-'*60}")
        print(f"  {'BIFO P@10':25s}  {s['bifo_p10']['mean']:>6.3f}  {s['bifo_p10']['sd']:>6.3f}  "
              f"{s['bifo_p10']['p25']:>5.3f}  {s['bifo_p10']['median']:>6.3f}  {s['bifo_p10']['p75']:>5.3f}")
        print(f"  {'BIFO Avg Precision':25s}  {ap['mean']:>6.3f}  {ap['sd']:>6.3f}  "
              f"{ap['p25']:>5.3f}  {ap['median']:>6.3f}  {ap['p75']:>5.3f}")
        print(f"  {'Rank improvement':25s}  {ri['mean']:>6.1f}  {ri['sd']:>6.1f}  "
              f"{ri['p25']:>5.1f}  {ri['median']:>6.1f}  {ri['p75']:>5.1f}")
        print(f"  Rank imp > 0: {s['rank_imp_positive']}/{s['n']} "
              f"({100*s['rank_imp_positive']/s['n']:.1f}%)")

    if primary_row:
        print(f"\n  Primary run ({len(primary_cuis)} seeds):")
        for k in ['bifo_p10','bifo_ap','rank_improvement']:
            print(f"    {k}: {primary_row[k]:.4f}")

    print(f"\n  Overall rank_imp > 0: {summary['overall']['rank_imp_positive']}/{len(boot_rows)} "
          f"({summary['overall']['rank_imp_positive_pct']:.1f}%)")
    print(f"\n  Output: {args.out_csv}")
    print(f"          {args.out_json}")


def main():
    p = argparse.ArgumentParser(
        description="Bootstrap resampling for KF cohort BIFO analyses"
    )
    p.add_argument('--kept-edges',   required=True)
    p.add_argument('--edges-merged', required=True)
    p.add_argument('--node-index',   required=True)
    p.add_argument('--bifo-scores',  default=None,
                   help='pathway_scores_standard.csv to freeze pathway universe')
    p.add_argument('--seed-pool',    required=True,
                   help='kf_chd_seeds.txt or kf_nbl_seeds.txt (gene symbols + counts)')
    p.add_argument('--ref-pathways', required=True,
                   help='Reference pathway CUI list (kf_chd_cilia_reference.txt)')
    p.add_argument('--out-csv',      required=True)
    p.add_argument('--out-json',     required=True)
    p.add_argument('--seed-sizes',   type=int, nargs='+', default=[10, 20, 30],
                   help='Seed set sizes to test (default: 10 20 30)')
    p.add_argument('--n-boots',      type=int, default=500,
                   help='Bootstrap draws per seed size (default: 500)')
    p.add_argument('--n-cores',      type=int, default=None)
    p.add_argument('--alpha',        type=float, default=0.5)
    p.add_argument('--min-members',  type=int, default=8)
    p.add_argument('--max-members',  type=int, default=300)
    p.add_argument('--random-seed',  type=int, default=42)
    run(p.parse_args())


if __name__ == '__main__':
    main()
