#!/usr/bin/env python3
"""Build a cilia/ciliopathy pathway reference from scored pathway universe."""
import argparse, csv, re
from pathlib import Path

CILIA_PATTERNS = [
    r'ciliopathies', r'ciliopathy', r'cilium', r'\bcilia\b', r'ciliary',
    r'ciliogenesis', r'intraflagellar', r'flagell', r'basal.body',
    r'hedgehog.signal', r'signaling.by.hedgehog', r'hedgehog.ligand',
    r'hedgehog.pathway', r'smoothened', r'left.right.asymmetry',
    r'nodal.flow', r'laterality', r'bardet.biedl', r'joubert',
    r'meckel', r'nephronophthisis', r'alstrom', r'oral.facial.digital',
    r'polycystic.kidney', r'primary.ciliary.dyskinesia',
]

def is_cilia(name):
    return any(re.search(p, name.lower()) for p in CILIA_PATTERNS)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scores',  required=True)
    ap.add_argument('--out',     required=True)
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    found = []
    with open(args.scores) as f:
        for row in csv.DictReader(f):
            name = row.get('name', '')
            cui  = row.get('concept_id', '')
            sab  = row.get('sab', '')
            if cui and is_cilia(name):
                found.append((cui, sab, name))

    found.sort(key=lambda x: (x[1], x[2]))
    with open(args.out, 'w') as f:
        f.write(f'# Cilia/ciliopathy pathway reference\n# {len(found)} pathways\n')
        for cui, sab, name in found:
            f.write(f'{cui}  # [{sab}] {name}\n')

    if args.verbose:
        print(f"Written: {args.out} ({len(found)} pathways)")
        for cui, sab, name in found:
            print(f"  {cui}  [{sab}] {name}")

if __name__ == '__main__':
    main()
