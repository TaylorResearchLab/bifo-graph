#!/usr/bin/env python3
"""
run_sm5_calibration_sweep.py
§SM5.3 — Calibration filter sensitivity sweep (CORRECTED VERSION)

Re-applies the `signal_to_null_mean` cutoff at multiple values to the
existing canonical pathway_scores_standard.csv files.

CANONICAL CALIBRATION RULE (from inspection of canonical CSV):
    A pathway is "calibrated" iff:
        signal_to_null_mean ≤ cutoff  (default cutoff = 10)
        AND null_z is non-NaN
    The boolean column `null_calibrated` in the CSV reflects this with
    the canonical cutoff = 10.

Pathways with signal_to_null_mean > cutoff are DEGENERATE: their null
distribution sits near zero and the null_z statistic is uninformative.
These are excluded from the calibrated set.

Sweep range: signal_to_null_mean cutoff ∈ {5, 7, 10, 15, 20}
Canonical cutoff: 10
Goal: quantify boundary instability discovered Apr 26 (~5% of pathways
flip across this cutoff between re-runs).

Verdict criterion (pre-specified):
    rank-within-calibrated Δ ≤ 5 across cutoffs ∈ [7, 15]
    where rank-within-calibrated is computed by sorting calibrated
    pathways by null_z descending.

Outputs:
    results/sm5_aggregate/calibration_sweep_kf_chd.json
    results/sm5_aggregate/calibration_sweep_kf_nbl.json

Usage:
    python3 pipeline/run_sm5_calibration_sweep.py
"""

import csv
import json
import sys
from pathlib import Path

# Pre-specified §SM5.3 design
CALIBRATION_CUTOFFS = [5, 7, 10, 15, 20]
CANONICAL_CUTOFF = 10
COHORTS = ['kf_chd', 'kf_nbl']

# Pathways tracked across the cilia cluster
TRACKED_PATHWAYS = [
    'WP_CILIOPATHIES',
    'WP_JOUBERT_SYNDROME',
    'WP_BARDET_BIEDL_SYNDROME',
    'REACTOME_ANCHORING_OF_THE_BASAL_BODY_TO_THE_PLASMA_MEMBRANE',
]

# Verdict thresholds (from agreed §SM5 criteria)
CALIBRATION_RANK_DELTA_MAX = 5
CALIBRATION_CUTOFF_RANGE = (7, 15)

# Paths
REPO_DIR = Path('/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph')
OUT_DIR = REPO_DIR / 'results' / 'sm5_aggregate'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _safe_float(x):
    if x is None or x == '':
        return None
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def load_pathway_scores(cohort):
    """Load pathway_scores_standard.csv for a cohort."""
    path = REPO_DIR / 'results' / cohort / 'pathway_scores_standard.csv'
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r['null_z_num'] = _safe_float(r.get('null_z'))
        r['signal_to_null_mean_num'] = _safe_float(r.get('signal_to_null_mean'))
        r['degree_norm_num'] = _safe_float(r.get('degree_norm'))
    return rows


def filter_calibrated(rows, cutoff):
    """Apply the calibration filter at a given cutoff.

    A pathway is calibrated iff:
      - signal_to_null_mean is non-NaN AND ≤ cutoff
      - null_z is non-NaN

    Then sort calibrated pathways by null_z DESCENDING (highest first)
    to assign rank-within-calibrated (rank 1 = highest null_z).
    """
    out = []
    for r in rows:
        s = r['signal_to_null_mean_num']
        z = r['null_z_num']
        if s is None or z is None:
            continue
        if s <= cutoff:  # ≤ is the canonical rule (not >)
            out.append(r)
    out.sort(key=lambda r: -r['null_z_num'])
    for i, r in enumerate(out, start=1):
        r['rank_within_calibrated'] = i
    return out


def find_pathway(filtered, name):
    for r in filtered:
        if r.get('name', '') == name:  # column is 'name', not 'pathway_name'
            return r
    return None


def boundary_flippers(rows, cutoff_low, cutoff_high):
    """Pathways with signal_to_null_mean ∈ (cutoff_low, cutoff_high]
    that flipped category between these two cutoffs.

    With ≤ rule: pathways IN at cutoff_high but OUT at cutoff_low.
    Specifically: signal_to_null_mean ∈ (cutoff_low, cutoff_high]
    (calibrated at cutoff_high but not at cutoff_low).
    """
    flippers = []
    for r in rows:
        s = r['signal_to_null_mean_num']
        if s is None:
            continue
        if cutoff_low < s <= cutoff_high:
            flippers.append({
                'name': r.get('name', ''),
                'signal_to_null_mean': s,
                'null_z': r['null_z_num'],
            })
    return flippers


def sweep_cohort(cohort):
    rows = load_pathway_scores(cohort)
    print(f'\n=== {cohort} ===')
    print(f'  Loaded {len(rows)} pathways from pathway_scores_standard.csv')

    out = {
        'cohort': cohort,
        'n_pathways_total': len(rows),
        'cutoffs_swept': CALIBRATION_CUTOFFS,
        'canonical_cutoff': CANONICAL_CUTOFF,
        'tracked_pathways': TRACKED_PATHWAYS,
        'calibration_rule': 'signal_to_null_mean <= cutoff AND null_z is not NaN',
        'rank_assignment_rule': 'sort calibrated set by null_z desc; rank 1 = highest null_z',
        'per_cutoff': [],
        'boundary_flippers_between_cutoffs': [],
    }

    for cutoff in CALIBRATION_CUTOFFS:
        filtered = filter_calibrated(rows, cutoff)
        cutoff_record = {
            'cutoff': cutoff,
            'n_calibrated': len(filtered),
            'tracked_ranks': {},
            'tracked_null_z': {},
        }
        for name in TRACKED_PATHWAYS:
            hit = find_pathway(filtered, name)
            if hit:
                cutoff_record['tracked_ranks'][name] = hit['rank_within_calibrated']
                cutoff_record['tracked_null_z'][name] = hit['null_z_num']
            else:
                cutoff_record['tracked_ranks'][name] = None
                cutoff_record['tracked_null_z'][name] = None
        out['per_cutoff'].append(cutoff_record)
        wp_rank = cutoff_record['tracked_ranks'].get('WP_CILIOPATHIES')
        wp_z = cutoff_record['tracked_null_z'].get('WP_CILIOPATHIES')
        z_str = f'{wp_z:.3f}' if wp_z is not None else 'NA'
        print(f'  cutoff={cutoff:>3}: n_calibrated={len(filtered):>5}  '
              f'WP_CILIOPATHIES rank-within-calibrated={wp_rank}  null_z={z_str}')

    # Boundary flippers
    for i in range(len(CALIBRATION_CUTOFFS) - 1):
        lo = CALIBRATION_CUTOFFS[i]
        hi = CALIBRATION_CUTOFFS[i + 1]
        flippers = boundary_flippers(rows, lo, hi)
        out['boundary_flippers_between_cutoffs'].append({
            'cutoff_low': lo,
            'cutoff_high': hi,
            'n_flippers': len(flippers),
            'flippers': flippers,
        })
        print(f'  boundary ({lo:>2},{hi:>2}]: '
              f'{len(flippers):>3} pathways flip between these cutoffs')

    # Verdict
    in_range = [
        c for c in out['per_cutoff']
        if CALIBRATION_CUTOFF_RANGE[0] <= c['cutoff'] <= CALIBRATION_CUTOFF_RANGE[1]
    ]
    cilio_ranks = [
        c['tracked_ranks'].get('WP_CILIOPATHIES')
        for c in in_range
        if c['tracked_ranks'].get('WP_CILIOPATHIES') is not None
    ]
    if cilio_ranks:
        rank_delta = max(cilio_ranks) - min(cilio_ranks)
        verdict = 'stable' if rank_delta <= CALIBRATION_RANK_DELTA_MAX else 'sensitive'
    else:
        rank_delta = None
        verdict = 'incomplete'

    out['verdict'] = {
        'criterion': f'rank-within-calibrated Δ ≤ {CALIBRATION_RANK_DELTA_MAX} '
                     f'across cutoffs ∈ [{CALIBRATION_CUTOFF_RANGE[0]}, '
                     f'{CALIBRATION_CUTOFF_RANGE[1]}]',
        'wp_ciliopathies_rank_delta': rank_delta,
        'wp_ciliopathies_ranks_in_range': cilio_ranks,
        'verdict': verdict,
    }
    print(f'  VERDICT (WP_CILIOPATHIES, cutoffs [{CALIBRATION_CUTOFF_RANGE[0]},'
          f'{CALIBRATION_CUTOFF_RANGE[1]}]): rank Δ = {rank_delta} → {verdict}')

    return out


def main():
    print('='*60)
    print('§SM5.3 calibration filter sensitivity sweep')
    print('='*60)
    print(f'Pre-specified cutoffs: {CALIBRATION_CUTOFFS}')
    print(f'Canonical cutoff: {CANONICAL_CUTOFF}')
    print(f'Calibration rule: signal_to_null_mean <= cutoff AND null_z is not NaN')
    print(f'Verdict criterion: rank Δ ≤ {CALIBRATION_RANK_DELTA_MAX} '
          f'across cutoffs ∈ [{CALIBRATION_CUTOFF_RANGE[0]}, '
          f'{CALIBRATION_CUTOFF_RANGE[1]}]')

    for cohort in COHORTS:
        result = sweep_cohort(cohort)
        out_path = OUT_DIR / f'calibration_sweep_{cohort}.json'
        with open(out_path, 'w') as f:
            json.dump(result, f, indent=2)
        print(f'  → wrote {out_path}')

    print('\nDone.')


if __name__ == '__main__':
    main()
