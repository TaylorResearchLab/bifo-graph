#!/usr/bin/env python3
"""
build_ncc_membership_edges.py

Converts NCC/cilia gene symbol pathway files into CUI-based pathway membership
edges compatible with the BIFO pipeline (same format as pathway_membership_edges.csv).

Also creates synthetic pathway concept nodes for each NCC/cilia set so they
can be scored alongside MSigDB/WikiPathways pathways in score_pathways.py.

Usage
-----
  python build_ncc_membership_edges.py \\
    --ncc-dir     ncc_cilia_pathways/ \\
    --nodes-csv   nodes.csv \\
    --out-edges   ncc_membership_edges.csv \\
    --out-nodes   ncc_pathway_nodes.csv \\
    --verbose

Output
------
  ncc_membership_edges.csv:
    source,target,predicate
    NCC:CARDIAC_NCC,C1234567,pathway_associated_with_gene
    ...
    (same format as pathway_membership_edges.csv — can be merged with it)

  ncc_pathway_nodes.csv:
    CUI,name,SAB
    NCC:CARDIAC_NCC,Cardiac neural crest subcircuit,NCC_CUSTOM
    ...
    (merge with nodes.csv so score_pathways.py can see pathway nodes)

Design notes
------------
  - NCC/cilia pathway nodes get synthetic CUIs in the "NCC:" namespace
    (e.g. NCC:CARDIAC_NCC) to avoid collision with UMLS CUIs
  - Gene membership uses predicate pathway_associated_with_gene 
    (same as MSigDB membership) so existing pipeline handles it correctly
  - Genes not found in nodes.csv are skipped with a warning
  - The SAB for NCC pathway nodes is "NCC_CUSTOM" so score_pathways.py
    can optionally filter by SAB
"""

import argparse

def removesuffix(s, suffix):
    return s[:-len(suffix)] if suffix and s.endswith(suffix) else s

import csv
import sys
from pathlib import Path
from collections import defaultdict


# ─── HGNC lookup (from nodes.csv) ─────────────────────────────────────────────

def build_hgnc_lookup(nodes_csv: str):
    """Return {gene_symbol: CUI} for gene nodes.

    Handles two formats from DDKG exports:
      - SAB='HGNC', name='SYMBOL gene'   (standard UBKG format)
      - SAB='HGNC', CODE='SYMBOL'        (code-based format)
      - name='SYMBOL gene' regardless of SAB (fallback)
    """
    lookup = {}
    n_total = 0
    with open(nodes_csv, encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        raw_rows = list(reader)

    if not raw_rows:
        print("  WARNING: nodes.csv appears empty")
        return lookup

    # Strip leading/trailing whitespace from all column names and values
    rows = [{k.strip(): v.strip() for k, v in row.items()} for row in raw_rows]
    sample_keys = list(rows[0].keys()) if rows else []

    for row in rows:
        n_total += 1
        # Support both old (CUI/SAB/CODE) and new (node_id/label/name/sab) formats
        sab  = (row.get('sab') or row.get('SAB') or '').strip()
        cui  = (row.get('node_id') or row.get('CUI') or '').strip()
        cui  = removesuffix(cui, ' CUI')
        name = row.get('name', '').strip()
        code = row.get('CODE', '').strip()

        # Only register HGNC nodes as gene symbols
        if sab != 'HGNC':
            continue

        # Method 1: name ends with " gene" — primary DDKG format
        if name.endswith(' gene'):
            sym = name[:-5].strip()
            if sym:
                lookup[sym] = cui
                lookup[sym.upper()] = cui

        # Method 2: CODE is plain gene symbol (fallback for old format)
        elif code and code != cui and not code.isdigit() and ' ' not in code:
            lookup[code] = cui
            lookup[code.upper()] = cui

    print(f"  {n_total} total nodes in nodes.csv")
    print(f"  {len(lookup)//2} gene symbols resolved")
    if sample_keys:
        print(f"  Columns: {sample_keys[:6]}")
    sample = [k for k in list(lookup.keys())[:40] if k == k.upper() and len(k) < 15]
    if sample:
        print(f"  Sample resolved symbols: {sample[:8]}")
    return lookup


# ─── Read NCC gene set files ───────────────────────────────────────────────────

def read_ncc_sets(ncc_dir: str):
    """
    Read all .txt files in ncc_dir.
    Returns {set_name: {'description': str, 'genes': [str]}}
    """
    sets = {}
    for fpath in sorted(Path(ncc_dir).glob('*.txt')):
        if fpath.name == 'MANIFEST.txt':
            continue
        name = fpath.stem
        desc = ''
        genes = []
        with open(fpath) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('# ') and not desc:
                    pass  # skip first comment (name)
                elif line.startswith('# ') and not genes:
                    candidate = line[2:].strip()
                    # second comment line is description
                    if not candidate.startswith('Source:') and not candidate.startswith('N_genes:') \
                       and not candidate.startswith('Format:'):
                        desc = candidate
                elif not line.startswith('#'):
                    genes.append(line)
        sets[name] = {'description': desc or name, 'genes': genes}
    return sets


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Build NCC/cilia pathway membership edges for BIFO pipeline"
    )
    p.add_argument('--ncc-dir',    required=True,
                   help='Directory containing NCC/cilia .txt gene set files')
    p.add_argument('--nodes-csv',  required=True,
                   help='nodes.csv from DDKG export (for HGNC CUI lookup)')
    p.add_argument('--out-edges',  required=True,
                   help='Output membership edges CSV')
    p.add_argument('--out-nodes',  required=True,
                   help='Output pathway concept nodes CSV (merge with nodes.csv)')
    p.add_argument('--verbose',    action='store_true')
    args = p.parse_args()

    print(f"Loading HGNC lookup from {args.nodes_csv}...")
    hgnc = build_hgnc_lookup(args.nodes_csv)
    print(f"  {len(hgnc)//2} HGNC genes in graph")

    print(f"Reading NCC gene sets from {args.ncc_dir}...")
    ncc_sets = read_ncc_sets(args.ncc_dir)
    print(f"  {len(ncc_sets)} gene sets")
    if len(ncc_sets) == 0:
        print(f"  WARNING: No .txt files found in {args.ncc_dir}")
        print(f"  Make sure the path points to the directory containing NPB_TF_CORE.txt etc.")
        print(f"  Absolute path tried: {Path(args.ncc_dir).resolve()}")

    # Build edges and collect stats
    edge_rows = []
    node_rows = []
    stats = {}

    for set_name, info in ncc_sets.items():
        pathway_cui = f"NCC:{set_name}"
        desc = info['description']
        genes = info['genes']

        found_cuis = []
        missing = []
        for sym in genes:
            cui = hgnc.get(sym) or hgnc.get(sym.upper())
            if cui:
                found_cuis.append((sym, cui))
                # pathway → gene edge
                edge_rows.append({
                    'source': pathway_cui,
                    'target': cui,
                    'predicate': 'pathway_associated_with_gene',
                })
                # gene → pathway edge (inverse)
                edge_rows.append({
                    'source': cui,
                    'target': pathway_cui,
                    'predicate': 'inverse_pathway_associated_with_gene',
                })
            else:
                missing.append(sym)

        # Pathway concept node
        node_rows.append({
            'node_id': pathway_cui,
            'label':   'Concept',
            'name':    desc,
            'sab':     'NCC_CUSTOM',
        })

        stats[set_name] = {
            'n_genes_in_set': len(genes),
            'n_found_in_graph': len(found_cuis),
            'n_missing': len(missing),
            'missing': missing,
        }

        if args.verbose or missing:
            pct = 100 * len(found_cuis) / len(genes) if genes else 0
            print(f"  {set_name:35s}  {len(found_cuis):3d}/{len(genes):3d} ({pct:.0f}%)"
                  + (f"  MISSING: {missing}" if missing else ""))

    # Write edges
    out_edges = Path(args.out_edges)
    with open(out_edges, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['source', 'target', 'predicate'])
        writer.writeheader()
        writer.writerows(edge_rows)
    print(f"\nWritten: {out_edges}  ({len(edge_rows)} edges, "
          f"{len(edge_rows)//2} gene-pathway pairs)")

    # Write nodes
    out_nodes = Path(args.out_nodes)
    with open(out_nodes, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['node_id', 'label', 'name', 'sab'])
        writer.writeheader()
        writer.writerows(node_rows)
    print(f"Written: {out_nodes}  ({len(node_rows)} pathway nodes)")

    # Summary
    print(f"\nCoverage summary:")
    total_gene_instances = sum(s['n_genes_in_set'] for s in stats.values())
    total_found = sum(s['n_found_in_graph'] for s in stats.values())
    all_missing = set(g for s in stats.values() for g in s['missing'])
    print(f"  Gene instances in sets:   {total_gene_instances}")
    pct = (100*total_found/total_gene_instances) if total_gene_instances > 0 else 0
    print(f"  Found in graph:           {total_found} ({pct:.0f}%)")
    print(f"  Unique genes missing:     {len(all_missing)}")
    if all_missing:
        print(f"  Missing symbols: {sorted(all_missing)[:20]}"
              + ("..." if len(all_missing) > 20 else ""))

    print(f"\nNext steps:")
    print(f"  1. Merge {out_nodes.name} into nodes.csv for your export:")
    print(f"     cat nodes.csv {out_nodes.name} > nodes_extended.csv")
    print(f"  2. Merge {out_edges.name} into your merged edges:")
    print(f"     (The run script will handle this if added to step [1/3])")
    print(f"  3. Run score_pathways.py with --chd-pathways pointing to")
    print(f"     a reference file listing NCC:* pathway CUIs you want evaluated")


if __name__ == '__main__':
    main()
