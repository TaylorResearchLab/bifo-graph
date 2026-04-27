#!/usr/bin/env python3
"""
compute_disease_relevant_intersections.py
§SM5.2 — Compute disease-relevant gene/CUI intersections for both cohorts.

Reads the disease-relevant reference files, maps gene symbols to UMLS CUIs
via the nodes_clean_noncc.csv files, intersects with the cohort seed CUI
files, and reports the disease-relevant pool sizes for §SM5.2.

Three pools computed:
  - KF-CHD (Jin 2017 H-CHD, 212 genes)
  - KF-NBL panel (Kim 2023 eTable 2 panel, 166 genes)
  - KF-NBL observed (Kim 2023 eTable 3 observed-variant subset, 54 genes)

Outputs intersected CUI lists to:
  data/cohorts/chd/kf_chd_disease_relevant_seed_cuis.txt
  data/cohorts/nbl/kf_nbl_disease_relevant_panel_seed_cuis.txt
  data/cohorts/nbl/kf_nbl_disease_relevant_observed_seed_cuis.txt

These intersected CUI files are the disease-relevant seed sets used in
§SM5.2 disease-relevant analysis. Pool sizes determine subset design:
  ≥200 → full 100×n=50/100/200 subset design works
  50-200 → smaller subset sizes (n=20-30) work
  <50 → single deterministic run on full intersection

Run from bifo-graph repo root:
    python3 pipeline/compute_disease_relevant_intersections.py
"""

import csv
import gzip
import sys
from pathlib import Path

REPO_DIR = Path('/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph')

POOLS = {
    'chd_jin2017': {
        'reference_file':    REPO_DIR / 'data/cohorts/chd/kf_chd_jin2017_reference.txt',
        'nodes_file':        REPO_DIR / 'results/kf_chd/nodes_clean_noncc.csv.gz',
        'seed_cuis_file':    REPO_DIR / 'data/cohorts/chd/kf_chd_seed_cuis.txt',
        'output_file':       REPO_DIR / 'data/cohorts/chd/kf_chd_disease_relevant_seed_cuis.txt',
        'cohort_label':      'KF-CHD',
        'pool_label':        'Jin 2017 H-CHD panel',
        'reference_source':  'Jin et al. 2017 Nat Genet, 212 H-CHD genes (Supp Data Set 2)',
    },
    'nbl_kim2023_panel': {
        'reference_file':    REPO_DIR / 'data/cohorts/nbl/kf_nbl_kim2023_panel_reference.txt',
        'nodes_file':        REPO_DIR / 'results/kf_nbl/nodes_clean_noncc.csv.gz',
        'seed_cuis_file':    REPO_DIR / 'data/cohorts/nbl/kf_nbl_seed_cuis.txt',
        'output_file':       REPO_DIR / 'data/cohorts/nbl/kf_nbl_disease_relevant_panel_seed_cuis.txt',
        'cohort_label':      'KF-NBL',
        'pool_label':        'Kim 2023 panel (eTable 2)',
        'reference_source':  'Kim et al. 2023 medRxiv, 166 CPG panel (eTable 2)',
    },
    'nbl_kim2023_observed': {
        'reference_file':    REPO_DIR / 'data/cohorts/nbl/kf_nbl_kim2023_observed_reference.txt',
        'nodes_file':        REPO_DIR / 'results/kf_nbl/nodes_clean_noncc.csv.gz',
        'seed_cuis_file':    REPO_DIR / 'data/cohorts/nbl/kf_nbl_seed_cuis.txt',
        'output_file':       REPO_DIR / 'data/cohorts/nbl/kf_nbl_disease_relevant_observed_seed_cuis.txt',
        'cohort_label':      'KF-NBL',
        'pool_label':        'Kim 2023 observed (eTable 3)',
        'reference_source':  'Kim et al. 2023 medRxiv, 54 observed-variant genes (eTable 3)',
    },
}


def load_reference_genes(path):
    """Load gene symbols from reference file. Skip comment lines."""
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            genes.add(line)
    return genes


def load_seed_cuis(path):
    """Load seed CUIs as a set. Skip comment lines (#) and blanks.

    Seed file format is one CUI per line, optionally followed by inline
    comments (e.g. "C1414995  # GATA4 [n_cases=1]"). We extract just the
    first whitespace-delimited token as the CUI.
    """
    cuis = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Take only the first token; ignore any inline comment
            cui = line.split()[0]
            cuis.add(cui)
    return cuis


# Cache for nodes file parsing — reused across pools sharing the same nodes file
_nodes_cache = {}


def build_gene_to_cui_map(nodes_path):
    """Parse nodes_clean_noncc.csv(.gz) and build a gene_symbol → CUI map.

    Cached per-path since both NBL pools share one nodes file.
    """
    nodes_path = str(nodes_path)
    if nodes_path in _nodes_cache:
        return _nodes_cache[nodes_path]

    open_fn = gzip.open if nodes_path.endswith('.gz') else open
    gene_to_cuis = {}

    with open_fn(nodes_path, 'rt') as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        print(f"  nodes file columns: {cols}", file=sys.stderr)

        cui_col = None
        for cand in ['cui', 'CUI', 'node_id', 'concept_id', 'id']:
            if cand in cols:
                cui_col = cand
                break

        sym_col = None
        for cand in ['gene_symbol', 'symbol', 'name', 'node_name', 'concept_name']:
            if cand in cols:
                sym_col = cand
                break

        sab_col = None
        for cand in ['sab', 'SAB', 'source']:
            if cand in cols:
                sab_col = cand
                break

        if cui_col is None or sym_col is None:
            raise RuntimeError(
                f"Could not find CUI and symbol columns in {nodes_path}. "
                f"Available: {cols}"
            )

        print(f"  Using CUI column: '{cui_col}', symbol column: '{sym_col}', "
              f"sab column: '{sab_col}'", file=sys.stderr)

        for row in reader:
            cui = row.get(cui_col, '').strip()
            sym = row.get(sym_col, '').strip()
            sab = row.get(sab_col, '').strip() if sab_col else None

            if not cui or not sym:
                continue

            # Strip " gene" suffix that DDKG sometimes adds
            if sym.endswith(' gene'):
                sym = sym[:-5]

            # Prefer HGNC-sourced rows; skip non-HGNC for gene resolution
            if sab_col and sab and sab != 'HGNC':
                continue

            gene_to_cuis.setdefault(sym, []).append(cui)

    _nodes_cache[nodes_path] = gene_to_cuis
    return gene_to_cuis


def intersect_pool(pool_key):
    """Run the intersection pipeline for one pool and return summary."""
    cfg = POOLS[pool_key]

    print(f"\n=== {cfg['cohort_label']} — {cfg['pool_label']} ===")
    print(f"  Reference: {cfg['reference_source']}")

    ref_genes = load_reference_genes(cfg['reference_file'])
    print(f"  Reference genes loaded: {len(ref_genes)}")

    seed_cuis = load_seed_cuis(cfg['seed_cuis_file'])
    print(f"  Cohort seed CUIs: {len(seed_cuis)}")

    print(f"  Building gene→CUI map from {cfg['nodes_file']}...")
    gene_to_cuis = build_gene_to_cui_map(cfg['nodes_file'])
    print(f"  Total gene symbols in nodes file: {len(gene_to_cuis)}")

    resolved = {}
    unresolved = []
    for gene in sorted(ref_genes):
        cuis = gene_to_cuis.get(gene)
        if cuis:
            resolved[gene] = cuis
        else:
            unresolved.append(gene)

    print(f"  Reference genes resolved to ≥1 CUI: {len(resolved)} / {len(ref_genes)}")
    if unresolved:
        print(f"  Unresolved reference genes ({len(unresolved)}):")
        for g in unresolved[:20]:
            print(f"    {g}")
        if len(unresolved) > 20:
            print(f"    ... and {len(unresolved) - 20} more")

    intersection = []
    for gene, cuis in resolved.items():
        for cui in cuis:
            if cui in seed_cuis:
                intersection.append((gene, cui))
                break

    print(f"  Reference genes in cohort seed list: {len(intersection)}")

    cfg['output_file'].parent.mkdir(parents=True, exist_ok=True)
    with open(cfg['output_file'], 'w') as f:
        f.write(f"# {cfg['cohort_label']} disease-relevant seed CUIs ({cfg['pool_label']})\n")
        f.write(f"# Source: {cfg['reference_source']}\n")
        f.write(f"# Intersection with {cfg['seed_cuis_file'].name}\n")
        f.write(f"# {len(intersection)} CUIs (from {len(ref_genes)} reference genes)\n")
        f.write(f"# Format: CUI\\tGENE_SYMBOL\n")
        for gene, cui in sorted(intersection):
            f.write(f"{cui}\t{gene}\n")
    print(f"  Wrote {cfg['output_file']}")

    return {
        'pool_key':       pool_key,
        'cohort':         cfg['cohort_label'],
        'pool':           cfg['pool_label'],
        'n_reference':    len(ref_genes),
        'n_resolved':     len(resolved),
        'n_unresolved':   len(unresolved),
        'n_intersection': len(intersection),
    }


def feasibility_verdict(n):
    """Translate intersection size to subset-design feasibility."""
    if n >= 200:
        return "full design (100×n=50/100/200)"
    elif n >= 100:
        return "subset design n=50/100 only"
    elif n >= 50:
        return "subset design n=20-30 only"
    elif n >= 10:
        return "single deterministic run only"
    else:
        return "POOL TOO SMALL (<10) — review"


def main():
    print("="*70)
    print("§SM5.2 disease-relevant gene intersection computation")
    print("="*70)

    summaries = []
    for pool_key in POOLS:
        try:
            summary = intersect_pool(pool_key)
            summaries.append(summary)
        except Exception as e:
            print(f"\n  ERROR for {pool_key}: {e}", file=sys.stderr)
            raise

    print("\n" + "="*70)
    print("Summary — disease-relevant seed pool sizes for §SM5.2 design")
    print("="*70)
    print(f"{'Pool':<35} {'Ref':>6} {'Resolv':>7} {'In_Seeds':>9} | Feasibility")
    print(f"{'':<35} {'':>6} {'':>7} {'':>9} |")
    for s in summaries:
        label = f"{s['cohort']:<8} {s['pool']}"
        verdict = feasibility_verdict(s['n_intersection'])
        print(f"{label:<35} {s['n_reference']:>6} {s['n_resolved']:>7} "
              f"{s['n_intersection']:>9} | {verdict}")

    print()
    print("Send this output to lock §SM5.2 subset design.")


if __name__ == '__main__':
    main()
