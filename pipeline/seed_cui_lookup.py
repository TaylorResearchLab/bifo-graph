#!/usr/bin/env python3
"""
seed_cui_lookup.py

Maps HGNC gene symbols to UMLS CUIs for use as BIFO pipeline seed nodes.

Two modes:
  1. Lookup against existing node index (uses nodes already in the graph export)
  2. Lookup against a provided HGNC→CUI mapping file (broader lookup)

Usage
-----
  # Mode 1: Look up against existing graph nodes.csv
  python seed_cui_lookup.py \\
    --gene-list   kf_chd_seeds.txt \\
    --nodes-csv   nodes.csv \\
    --out         kf_chd_seed_cuis.txt

  # Mode 2: Look up against HGNC mapping table (nodes.csv from any DDKG export)
  python seed_cui_lookup.py \\
    --gene-list   kf_chd_seeds.txt \\
    --nodes-csv   nodes.csv \\
    --out         kf_chd_seed_cuis.txt \\
    --verbose

Input gene-list format
----------------------
One gene symbol per line, optional count column (tab-separated):
  GATA4
  NKX2-5
  TBX5   15       ← n_cases column (used for ranking, not required)
  DNAH11  8
  # comments ignored

Output format
-------------
  C1414995  # GATA4  [n_cases=15]
  C1413784  # NKX2-5
  (suitable for direct use as --seed-nodes in bifo_conditioning.py)
"""

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

import csv
import sys
from pathlib import Path




def read_gene_list(path: str):
    """Read gene symbols (and optional counts) from input file."""
    genes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            symbol = parts[0].strip()
            count = int(parts[1]) if len(parts) > 1 and parts[1].strip().isdigit() else None
            if symbol:
                genes.append((symbol, count))
    return genes


def build_hgnc_lookup(nodes_csv: str):
    """
    Build {gene_symbol: CUI} from nodes.csv.
    
    In DDKG nodes.csv, HGNC gene concepts have:
      - SAB = 'HGNC'
      - name = 'SYMBOL gene' (e.g. 'GATA4 gene')
      OR CODE = HGNC symbol
    
    We try both formats.
    """
    lookup = {}  # symbol → CUI
    cui_to_name = {}
    
    with open(nodes_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Strip whitespace from all keys and values
            row = {k.strip(): v.strip() for k, v in row.items()}
            # Support both old (CUI/SAB/CODE) and new (node_id/sab/name) formats
            sab  = (row.get('sab') or row.get('SAB') or '').strip()
            cui  = (row.get('node_id') or row.get('CUI') or '').strip()
            cui  = removesuffix(cui, ' CUI')
            name = row.get('name', '').strip()
            code = row.get('CODE', '').strip()

            if sab != 'HGNC':
                continue

            # Method 1: name = "SYMBOL gene" (confirmed DDKG format)
            if name.endswith(' gene'):
                symbol = name[:-5].strip()
                lookup[symbol] = cui
                lookup[symbol.upper()] = cui

            # Method 2: CODE = HGNC symbol (old format fallback)
            if code and code != cui:
                lookup[code] = cui
                lookup[code.upper()] = cui

            cui_to_name[cui] = name

    return lookup, cui_to_name


def main():
    p = argparse.ArgumentParser(description="Map HGNC gene symbols to UMLS CUIs")
    p.add_argument('--gene-list', required=True,
                   help='Gene symbol list (one per line, optional count column)')
    p.add_argument('--nodes-csv', required=True,
                   help='nodes.csv from DDKG export (must include HGNC SAB entries)')
    p.add_argument('--out', required=True,
                   help='Output seed CUI file (suitable for --seed-nodes)')
    p.add_argument('--verbose', action='store_true',
                   help='Print lookup details')
    p.add_argument('--min-count', type=int, default=None,
                   help='Minimum case count to include (requires count column in gene-list)')
    p.add_argument('--top-n', type=int, default=None,
                   help='Select top N genes by count (requires count column)')
    args = p.parse_args()
    
    print(f"Reading gene list: {args.gene_list}")
    genes = read_gene_list(args.gene_list)
    print(f"  {len(genes)} genes loaded")
    
    if args.min_count is not None:
        before = len(genes)
        genes = [(s, c) for s, c in genes if c is None or c >= args.min_count]
        print(f"  min_count filter ({args.min_count}): {before} → {len(genes)}")
    
    if args.top_n is not None:
        genes_with_count = [(s, c or 0) for s, c in genes]
        genes_with_count.sort(key=lambda x: -x[1])
        genes = genes_with_count[:args.top_n]
        print(f"  top_n ({args.top_n}): selected top {len(genes)} by count")
    
    print(f"Building HGNC lookup from: {args.nodes_csv}")
    lookup, cui_to_name = build_hgnc_lookup(args.nodes_csv)
    print(f"  {len(lookup)//2} HGNC genes in graph")
    
    # Look up each gene
    found, missing = [], []
    for symbol, count in genes:
        cui = lookup.get(symbol) or lookup.get(symbol.upper())
        if cui:
            found.append((symbol, cui, count, cui_to_name.get(cui, '')))
        else:
            missing.append((symbol, count))
    
    print(f"\nResults: {len(found)} found, {len(missing)} not in graph")
    
    if missing:
        print(f"\nMissing from graph:")
        for sym, count in missing:
            count_str = f"  (n={count})" if count else ""
            print(f"  {sym}{count_str}")
    
    if args.verbose:
        print(f"\nFound:")
        for sym, cui, count, name in found:
            count_str = f"  n={count}" if count else ""
            print(f"  {sym:15s} → {cui:12s}  ({name}){count_str}")
    
    # Write output
    out_path = Path(args.out)
    with open(out_path, 'w') as f:
        f.write(f"# Cohort variant-derived seeds\n")
        f.write(f"# Source: {args.gene_list}\n")
        f.write(f"# Graph: {args.nodes_csv}\n")
        f.write(f"# Found: {len(found)}/{len(genes)} genes\n")
        if args.min_count:
            f.write(f"# min_count filter: {args.min_count}\n")
        if args.top_n:
            f.write(f"# top_n: {args.top_n}\n")
        f.write("#\n")
        for sym, cui, count, name in found:
            count_str = f"  # {sym}" + (f" [n_cases={count}]" if count else f"  # {sym}")
            f.write(f"{cui}{count_str}\n")
    
    print(f"\nWritten: {out_path}  ({len(found)} CUIs)")
    
    # Summary statistics
    if any(c for _, c in genes if c):
        counts = sorted([c for _, c in genes if c], reverse=True)
        print(f"\nCount distribution of selected genes:")
        print(f"  max={counts[0]}  median={counts[len(counts)//2]}  min={counts[-1]}")
        for threshold in [10, 5, 3, 2, 1]:
            n = sum(1 for c in counts if c >= threshold)
            print(f"  n_cases >= {threshold}: {n} genes")


if __name__ == '__main__':
    main()
