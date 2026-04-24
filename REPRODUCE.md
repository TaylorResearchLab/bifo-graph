# Reproducing All Analyses in Taylor et al. 2026

This document provides exact commands to reproduce every result reported in the manuscript. All analyses were run using the frozen pipeline in this repository. Pre-computed results are archived in `results/` and can be verified without re-running.

## Prerequisites

```bash
git clone https://github.com/TaylorResearchLab/bifo-graph
cd bifo-graph
pip install numpy pandas scipy scikit-learn pyyaml
```

Python 3.7+ required. All analyses use `pipeline/bifo_conditioning.py` and `pipeline/score_pathways.py` from this repository. The frozen configuration is `config/bifo_mapping.yaml` (v0.7.1).

---

## Analysis 1 — CHD Curated Benchmark (Tables 1–4, Figures 1–4)

### Step 1.1 — BIFO Conditioning and PPR Propagation

```bash
# Merge raw + membership edges
python3 -c "
import pandas as pd
raw = pd.read_csv('data/benchmark/chd_curated_edges_raw.csv.zip')
mem = pd.read_csv('data/benchmark/chd_curated_pathway_membership_edges.csv.zip')
pd.concat([raw, mem], ignore_index=True).to_csv('results/chd_benchmark/edges_merged_chd.csv', index=False)
print(f'Merged: {len(raw)} + {len(mem)} = {len(raw)+len(mem)} edges')
"

# Run conditioning (four-arm PPR)
python3 pipeline/bifo_conditioning.py \
  --nodes   data/benchmark/chd_curated_nodes.csv \
  --edges   results/chd_benchmark/edges_merged_chd.csv \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --out-json results/chd_benchmark/results_full.json
```

**Frozen outputs:** `results/chd_benchmark/results_full.json`, `results_full_scores_cond.npy`, `results_full_scores_raw.npy`, `results_full_kept_edges.csv.zip`, `results_node_index.json`

**Expected:** alpha=0.5, 105,192 kept edges (57,005 propagating), conditioned AUPRC=0.1947, entropy=4.934

### Step 1.2 — Ablation Arm (edges_raw only, no bridge edges)

```bash
python3 pipeline/bifo_conditioning.py \
  --nodes   data/benchmark/chd_curated_nodes.csv \
  --edges   data/benchmark/chd_curated_edges_raw.csv.zip \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --out-json results/chd_benchmark/results_ablation.json
```

### Step 1.3 — Mechanistic-Only Arm

```bash
python3 pipeline/bifo_conditioning.py \
  --nodes   data/benchmark/chd_curated_nodes.csv \
  --edges   results/chd_benchmark/edges_merged_chd.csv \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --mechanistic-only \
  --out-json results/chd_benchmark/results_mech.json
```

### Step 1.4 — Pathway Scoring: Full BIFO Arm (Table 3, row 1)

```bash
# Extract kept_edges from zip first
python3 -c "
import zipfile, pandas as pd
z = zipfile.ZipFile('results/chd_benchmark/results_full_kept_edges.csv.zip')
name = [n for n in z.namelist() if '/' not in n][0]
pd.read_csv(z.open(name), low_memory=False).to_csv('results/chd_benchmark/chd_kept_edges.csv', index=False)
"

python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         results/chd_benchmark/edges_merged_chd.csv \
  --edges-conditioned results/chd_benchmark/chd_kept_edges.csv \
  --scores-cond       results/chd_benchmark/results_full_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_full_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --out-csv  results/chd_benchmark/pathway_scores_full.csv \
  --out-json results/chd_benchmark/pathway_metrics_full.json
```

**Expected:** 550 pathways, P@10=0.70, enrichment=21.4×, mean rank=86.6, rank improvement=+125.4

### Step 1.5 — Pathway Scoring with Empirical Null (Methods §5)

The null model is run as part of Step 1.4 by adding `--n-permutations 1000 --null-type membership-rewiring`. The null results are written to `pathway_scores_full.csv` alongside the scoring results. To reproduce the null results in isolation:

```bash
python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         results/chd_benchmark/edges_merged_chd.csv \
  --edges-conditioned results/chd_benchmark/chd_kept_edges.csv \
  --scores-cond       results/chd_benchmark/results_full_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_full_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 120 \
  --out-csv  results/chd_benchmark/pathway_scores_full.csv \
  --out-json results/chd_benchmark/pathway_metrics_full.json
```

**Expected (pathway-node rewiring null):** 16/550 pathways significant at q<0.05; BRUNEAU q=0.034 null_z=15.6

**Expected (member-level stratified null):** 17/550 pathways significant at q<0.05; BRUNEAU null_z=10.7 q=0.032; WP_HEART_DEVELOPMENT null_z=10.1 q=0.032

### Step 1.6 — Pathway Scoring: Ablation and Mechanistic Arms (Table 3)

```bash
# Ablation arm
python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         results/chd_benchmark/edges_merged_chd.csv \
  --edges-conditioned results/chd_benchmark/results_ablation_kept_edges.csv.zip \
  --scores-cond       results/chd_benchmark/results_ablation_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_ablation_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --out-csv  results/chd_benchmark/pathway_scores_ablation.csv \
  --out-json results/chd_benchmark/pathway_metrics_ablation.json

# Mechanistic-only arm
python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         results/chd_benchmark/edges_merged_chd.csv \
  --edges-conditioned results/chd_benchmark/results_mech_kept_edges.csv.zip \
  --scores-cond       results/chd_benchmark/results_mech_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_mech_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --out-csv  results/chd_benchmark/pathway_scores_mech.csv \
  --out-json results/chd_benchmark/pathway_metrics_mech.json
```

### Step 1.7 — Baseline Enrichment Comparison (Table 4)

```bash
python3 pipeline/baseline_enrichment.py \
  --edges-merged  results/chd_benchmark/edges_merged_chd.csv \
  --node-index    results/chd_benchmark/results_node_index.json \
  --scores-raw    results/chd_benchmark/results_full_scores_raw.npy \
  --scores-cond   results/chd_benchmark/results_full_scores_cond.npy \
  --bifo-scores   results/chd_benchmark/pathway_scores_full.csv \
  --chd-pathways  data/benchmark/chd_pathway_reference.txt \
  --seed-nodes    data/benchmark/chd_seed_nodes.txt \
  --out-csv       results/chd_benchmark/baseline_comparison.csv \
  --out-json      results/chd_benchmark/baseline_comparison.json
```

### Step 1.8 — C4 Pathway-Split Controls (Table 5)

```bash
# First extract kept_edges from zips
for ARM in c4_notch c4_mapk; do
  python3 -c "
import zipfile, pandas as pd, sys
arm = sys.argv[1]
z = zipfile.ZipFile(f'results/chd_benchmark/{arm}/results_kept_edges.csv.zip')
name = [n for n in z.namelist() if '/' not in n][0]
pd.read_csv(z.open(name), low_memory=False).to_csv(f'results/chd_benchmark/{arm}_kept_edges.csv', index=False)
" $ARM
done

# Score each C4 arm
for ARM in c4_notch c4_mapk; do
  python3 pipeline/score_pathways.py \
    --nodes             data/benchmark/chd_curated_nodes.csv \
    --edges-raw         results/chd_benchmark/edges_merged_chd.csv \
    --edges-conditioned results/chd_benchmark/${ARM}_kept_edges.csv \
    --scores-cond       results/chd_benchmark/${ARM}/results_scores_cond.npy \
    --scores-raw        results/chd_benchmark/${ARM}/results_scores_raw.npy \
    --node-index        results/chd_benchmark/${ARM}/results_node_index.json \
    --seed-nodes        data/benchmark/${ARM}_seed_nodes.txt \
    --chd-pathways      data/benchmark/${ARM}_pathway_reference.txt \
    --min-members 8 --max-members 300 \
    --out-csv  results/chd_benchmark/${ARM}/pathway_scores.csv \
    --out-json results/chd_benchmark/${ARM}/pathway_metrics.json
done
```

### Step 1.9 — Exhaustive Resampling (Table 6, Figure 6)

```bash
python3 pipeline/chd_resampling_exhaustive.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-merged      results/chd_benchmark/edges_merged_chd.csv \
  --mapping           config/bifo_mapping.yaml \
  --gene-pool         data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --out-csv           results/chd_benchmark/resampling_results.csv \
  --out-json          results/chd_benchmark/resampling_summary.json \
  --n-cores           120
```

**Expected:** 3,003 splits, 100% positive rank improvement, 94.4% P@10≥0.30

---

## Analysis 2 — KF-CHD Cohort (Results §8, Methods §10–15)

### Data requirements

KF-CHD graph data files (`edges_all_noncc.csv.gz`, `nodes_clean_noncc.csv.gz`) are archived in `results/kf_chd/`. Seed CUIs are in `data/cohorts/chd/kf_chd_seed_cuis.txt` (1,276 CUIs; 1,276/1,287 resolved).

```bash
# Already exists: kf_chd_edges_merged.csv
# Produced by merging kf_chd_edges_raw_clean.csv + kf_chd_pathway_membership_edges_clean.csv
```

Seed file: `kf_chd_seed_cuis.txt` (1,276 CUIs; derived from `data/cohorts/chd/kf_chd_seeds_maf001.txt`, 1,276/1,287 resolved to graph CUIs)

### Step 2.1 — Conditioning and Propagation

```bash


python3 pipeline/bifo_conditioning.py \
  --nodes   results/kf_chd/nodes_clean_noncc.csv.gz \
  --edges   results/kf_chd/edges_all_noncc.csv.gz \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/cohorts/chd/kf_chd_seed_cuis.txt \
  --heldout-nodes data/cohorts/chd/kf_chd_seed_cuis.txt \
  --out-json results/kf_chd/results.json
```

**Frozen outputs:** `results/kf_chd/`

**Expected:** 2,490,510 kept edges (44.1% of concept edges), 1,276/1,276 seeds resolved. Bridge edges: 1,026,968 / 2,482,752 propagating = 41.4%

### Step 2.2 — Pathway Scoring (Standard Universe)

```bash
python3 pipeline/score_pathways.py \
  --nodes             results/kf_chd/nodes_clean_noncc.csv.gz \
  --edges-raw         results/kf_chd/edges_all_noncc.csv.gz \
  --edges-conditioned results/kf_chd/results_kept_edges.csv.zip \
  --scores-cond       results/kf_chd/results_scores_cond.npy \
  --scores-raw        results/kf_chd/results_scores_raw.npy \
  --node-index        results/kf_chd/results_node_index.json \
  --seed-nodes        data/cohorts/chd/kf_chd_seed_cuis.txt \
  --chd-pathways      data/cohorts/chd/kf_chd_cilia_reference.txt \
  --allowed-name-prefixes HALLMARK_ REACTOME_ WP_ KEGG_ BIOCARTA_ PID_ \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 128 \
  --out-csv  results/kf_chd/pathway_scores_standard.csv \
  --out-json results/kf_chd/pathway_metrics_standard.json
```

**Expected (scoring):** WP_CILIOPATHIES rank 43/2,130 (top 2%), degree_norm=8.509e-06; WP_JOUBERT_SYNDROME rank 34 (null degenerate in CHD). Pathway universe: 2,130 (CGP sets excluded via --allowed-name-prefixes).

**Expected (pathway-node rewiring null):** WP_CILIOPATHIES null_z=41.19, empirical q=0.0081 (valid — bridge edges 41.4% of propagating graph; see Methods §5)

**Expected (member-level null):** WP_CILIOPATHIES member_mean null_z=1.39, empirical p=0.057 (not significant at member level; signal concentrated at pathway node)

> Note: `score_pathways.py` runs both nulls in a single invocation when `--null-type membership-rewiring` is specified. Output columns: `empirical_q`, `null_z` (pathway-node rewiring null); `member_mean_q`, `member_mean_null_z` (member-level null).

### Step 2.3 — Baseline Enrichment

```bash
python3 pipeline/baseline_enrichment.py \
  --edges-merged  results/kf_chd/edges_all_noncc.csv.gz \
  --node-index    results/kf_chd/results_node_index.json \
  --scores-raw    results/kf_chd/results_scores_raw.npy \
  --scores-cond   results/kf_chd/results_scores_cond.npy \
  --bifo-scores   results/kf_chd/pathway_scores_standard.csv \
  --chd-pathways  data/cohorts/chd/kf_chd_cilia_reference.txt \
  --seed-nodes    data/cohorts/chd/kf_chd_seed_cuis.txt \
  --small-universe \
  --out-csv  results/kf_chd/baseline_comparison.csv \
  --out-json results/kf_chd/baseline_comparison.json
```

**Expected:** WP_CILIOPATHIES rank 1 under seed Fisher (corrected, BH p=4.53×10⁻³¹, overlap=61/170); Fisher P@10=0.20, AP=0.130; BIFO P@10=0.00, AP=0.012, mean_ref_rank=889; raw PPR GSEA rank 1,994; cond PPR GSEA rank 456; neighborhood Fisher P@10=0.0

---

## Analysis 3 — KF-NBL Cohort (Results §8, Methods §10–15)

Independent neuroblastoma cohort (554 probands, dbGaP phs001436). Same pipeline as KF-CHD. Seed CUIs: `data/cohorts/nbl/kf_nbl_seed_cuis.txt` (1,395 CUIs; 1,395/1,406 resolved). Frozen outputs in `results/kf_nbl/`.

Seed file: `kf_nbl_seed_cuis.txt` (1,395 CUIs; NBL rare variant carriers at MAF ≤ 0.001)

### Step 3.1 — Conditioning and Propagation

```bash


python3 pipeline/bifo_conditioning.py \
  --nodes   results/kf_nbl/nodes_clean_noncc.csv.gz \
  --edges   results/kf_nbl/edges_all_noncc.csv.gz \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/cohorts/nbl/kf_nbl_seed_cuis.txt \
  --heldout-nodes data/cohorts/nbl/kf_nbl_seed_cuis.txt \
  --out-json results/kf_nbl/results.json
```

**Expected:** 2,654,867 kept edges (45.3% of concept edges), 1,395/1,406 seeds resolved

### Step 3.2 — Pathway Scoring

```bash
python3 pipeline/score_pathways.py \
  --nodes             results/kf_nbl/nodes_clean_noncc.csv.gz \
  --edges-raw         results/kf_nbl/edges_all_noncc.csv.gz \
  --edges-conditioned results/kf_nbl/results_kept_edges.csv.zip \
  --scores-cond       results/kf_nbl/results_scores_cond.npy \
  --scores-raw        results/kf_nbl/results_scores_raw.npy \
  --node-index        results/kf_nbl/results_node_index.json \
  --seed-nodes        data/cohorts/nbl/kf_nbl_seed_cuis.txt \
  --chd-pathways      data/cohorts/nbl/kf_nbl_cilia_reference.txt \
  --allowed-name-prefixes HALLMARK_ REACTOME_ WP_ KEGG_ BIOCARTA_ PID_ \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 128 \
  --out-csv  results/kf_nbl/pathway_scores_standard.csv \
  --out-json results/kf_nbl/pathway_metrics_standard.json
```

**Expected (scoring):** WP_CILIOPATHIES rank 3/2,196, degree_norm=4.241e-06, null_z=18.37, q=0.014; WP_JOUBERT_SYNDROME rank 17. Pathway universe: 2,196 (CGP sets excluded).

**Expected (pathway-node rewiring null):** WP_CILIOPATHIES null_z=18.37, empirical q=0.014 (valid)

**Expected (member-level null):** WP_CILIOPATHIES member_mean null_z=2.43, empirical p=0.003 (significant — signal distributed across member genes in NBL)

### Step 3.3 — Baseline Enrichment

```bash
python3 pipeline/baseline_enrichment.py \
  --edges-merged  results/kf_nbl/edges_all_noncc.csv.gz \
  --node-index    results/kf_nbl/results_node_index.json \
  --scores-raw    results/kf_nbl/results_scores_raw.npy \
  --scores-cond   results/kf_nbl/results_scores_cond.npy \
  --bifo-scores   results/kf_nbl/pathway_scores_standard.csv \
  --chd-pathways  data/cohorts/nbl/kf_nbl_cilia_reference.txt \
  --seed-nodes    data/cohorts/nbl/kf_nbl_seed_cuis.txt \
  --small-universe \
  --out-csv  results/kf_nbl/baseline_comparison.csv \
  --out-json results/kf_nbl/baseline_comparison.json
```

**Expected:** WP_CILIOPATHIES rank 1 under seed Fisher (corrected, BH p=3.85×10⁻³²); Fisher P@10=0.20, AP=0.099; BIFO P@10=0.10, AP=0.044, mean_ref_rank=550; raw PPR GSEA rank 1,707; cond PPR GSEA rank 1,543; convergent with KF-CHD

---

## Verifying Pre-Computed Results used in the BIFO proof of principle paper.

All frozen outputs are in `results/`. To verify all reported numbers match the manuscript,
run from the `bifo-graph/` repo root:

```bash

python3 -c "
import json, pandas as pd

print('=== TABLE 1: Conditioning ===')
r = json.load(open('results/chd_benchmark/results_full.json'))
cov = r['coverage']
gs = r['graph_stats']
print(f'  Nodes: {cov[\"total_nodes\"]} (expected 34,523)')
print(f'  Kept edges: {cov[\"kept_edges\"]} (expected 105,192)')
print(f'  Pathway Contribution (propagating gene->PW): {gs[\"flow_class_distribution\"][\"Pathway Contribution\"]} (expected 43,698 = 76.7% of 57,005 propagating edges)')
print(f'  Cond AUPRC: {r[\"conditioned\"][\"auprc\"]:.4f} (expected 0.1947; paper Table 1 and results_full.json are now from the same run)')
print(f'  Cond entropy: {r[\"conditioned\"][\"entropy\"]:.3f} (expected 4.934)')
print(f'  Raw AUPRC: {r[\"raw\"][\"auprc\"]:.4f} (expected 0.2215)')
print(f'  Raw entropy: {r[\"raw\"][\"entropy\"]:.3f} (expected 5.728)')

print()
print('=== TABLE 3: Three-arm ablation ===')
for arm, ep10, eri in [('full',0.70,+125.4),('ablation',0.60,+13.5),('mech',0.00,34.7)]:
    r = json.load(open(f'results/chd_benchmark/pathway_metrics_{arm}.json'))
    p10 = r['metrics']['top10_precision']
    ri = r['metrics']['rank_improvement_cond_vs_raw']
    status = '✅' if abs(p10 - ep10) < 0.001 else '❌'
    print(f'  {status} {arm}: P@10={p10:.2f} rank_imp={ri:.1f} (expected P@10={ep10} rank_imp={eri})')

print()
print('=== TABLE 4: Baseline comparison ===')
b = json.load(open('results/chd_benchmark/baseline_comparison.json'))
expected = {'degree_overlap':0.40,'seed_fisher':0.30,'neighborhood_fisher':0.00,
            'raw_ppr_gsea':0.10,'cond_ppr_gsea':0.10,'bifo_full':0.70}
ap_expected = {'degree_overlap':0.343,'seed_fisher':0.156,'neighborhood_fisher':0.037,
               'raw_ppr_gsea':0.117,'cond_ppr_gsea':0.088,'bifo_full':0.408}
for m in b['methods']:
    meth = m['method']
    p10 = m.get('precision_at_10', 0)
    ap = m.get('average_precision', 0)
    ep10 = expected.get(meth, '?')
    eap = ap_expected.get(meth, '?')
    status = '✅' if ep10 == '?' or abs(p10 - ep10) < 0.001 else '❌'
    print(f'  {status} {meth}: P@10={p10:.2f} AP={ap:.3f} (expected P@10={ep10} AP={eap})')

print()
print('=== TABLE 5: C4 Controls ===')
for arm, ep10, eap, eri in [('c4_notch',0.50,0.450,-31.5),('c4_mapk',0.10,0.174,-54.4)]:
    m = json.load(open(f'results/chd_benchmark/{arm}/pathway_metrics.json'))
    b = json.load(open(f'results/chd_benchmark/{arm}/baseline_comparison.json'))
    p10 = m['metrics']['top10_precision']
    ri  = m['metrics']['rank_improvement_cond_vs_raw']
    ap  = next(x['average_precision'] for x in b['methods'] if x['method'] == 'bifo_full')
    p10_ok = '✅' if abs(p10 - ep10) < 0.001 else '❌'
    ap_ok  = '✅' if abs(ap - eap) < 0.005 else '❌'
    ri_ok  = '✅' if abs(ri - eri) < 0.2 else '❌'
    print(f'  {p10_ok} {arm}: P@10={p10:.2f} {ap_ok} AP={ap:.3f} {ri_ok} rank_imp={ri:.1f}')
    print(f'         (expected P@10={ep10} AP={eap} rank_imp={eri})')

print()
print('=== TABLE 6: Resampling ===')
r = json.load(open('results/chd_benchmark/resampling_summary.json'))
rob = r['robustness']
mets = r['metrics']
print(f'  Splits: {r[\"n_splits_total\"]} (expected 3,003)')
print(f'  P@10>=0.30: {rob[\"bifo_p10_ge_0.3\"]}/3003 = {rob[\"bifo_p10_ge_0.3\"]/3003*100:.1f}% (expected 94.4%)')
print(f'  P@10>=0.50: {rob[\"bifo_p10_ge_0.5\"]}/3003 = {rob[\"bifo_p10_ge_0.5\"]/3003*100:.1f}% (expected 43.8%)')
print(f'  Positive rank imp: {rob[\"rank_imp_positive\"]}/3003 (expected 3003)')
print(f'  Median P@10: {mets[\"bifo_p10\"][\"median\"]:.2f} (expected 0.40)')
print(f'  AP range: {mets[\"bifo_ap\"][\"min\"]:.3f}-{mets[\"bifo_ap\"][\"max\"]:.3f} (expected 0.136-0.448)')

print()
print('=== NULL MODEL: Benchmark (from pathway_scores_full.csv) ===')
df = pd.read_csv('results/chd_benchmark/pathway_scores_full.csv')
df = df.sort_values('degree_norm', ascending=False)
sig = (df.empirical_q < 0.05).sum()
sig_mm = (df.member_mean_q < 0.05).sum()
bruneau = df[df.name == 'BRUNEAU_SEPTATION_VENTRICULAR'].iloc[0]
print(f'  Rewiring null sig q<0.05: {sig} (expected 16)')
print(f'  BRUNEAU rewiring null_z: {bruneau.null_z:.1f} (expected 15.6)')
print(f'  BRUNEAU q: {bruneau.empirical_q:.3f} (expected 0.034)')
print(f'  Member null sig q<0.05: {sig_mm} (expected 17)')
print(f'  BRUNEAU member null_z: {bruneau.member_mean_null_z:.1f} (expected 10.7)')
print(f'  BRUNEAU member q: {bruneau.member_mean_q:.3f} (expected 0.032)')

wph = df[df.name == 'WP_HEART_DEVELOPMENT'].iloc[0]
wpc = df[df.name == 'WP_CARDIAC_PROGENITOR_DIFFERENTIATION'].iloc[0]
print(f'  WP_HEART null_z: {wph.null_z:.1f} (expected 16.3)')
print(f'  WP_CARDIAC_PROGENITOR null_z: {wpc.null_z:.1f} (expected 12.3)')

print()
print('=== KF-CHD GSEA ranks (from results/kf_chd/baseline_comparison.csv) ===')
df4 = pd.read_csv('results/kf_chd/baseline_comparison.csv')
for method, expected_rank in [('raw_ppr_gsea',1994),('cond_ppr_gsea',456),('neighborhood_fisher',2126)]:
    sub = df4[df4.method == method]
    cilia4 = sub[sub.name.str.contains('CILIOPATHIES', na=False)]
    rank = int(cilia4.iloc[0]['rank']) if len(cilia4) > 0 else 'not found'
    status = '✅' if rank == expected_rank else '❌'
    print(f'  {status} {method}: WP_CILIOPATHIES rank={rank} (expected {expected_rank})')
"
```

**Note on Cond AUPRC:** `results_full.json` and paper Table 1 both report 0.1947. These are now from the same clean run of the unidirectional pipeline.

For KF-CHD and KF-NBL verification, also run these commands from inside the bifo-graph root directory: 

```bash

python3 -c "
import json, pandas as pd

print('=== KF-CHD RESULTS ===')
rs = json.load(open('results/kf_chd/resampling_summary.json'))
pr = rs['primary_run']
print(f'  Primary run P@10: {pr[\"bifo_p10\"]} (expected 0.00)')
print(f'  Primary run AP: {pr[\"bifo_ap\"]:.3f} (expected 0.012)')
print(f'  Fisher P@10: {pr[\"sf_p10\"]} (expected 0.20)')
print(f'  Fisher AP: {pr[\"sf_ap\"]:.3f} (expected 0.130)')

# pathway_scores_standard.csv has ranking; null results are in separate file
df = pd.read_csv('results/kf_chd/pathway_scores_standard.csv')
df = df.sort_values('degree_norm', ascending=False)
df.index = range(1, len(df)+1)
cilia_rank = df[df.name == 'WP_CILIOPATHIES'].index[0]
print(f'  WP_CILIOPATHIES rank: {cilia_rank} (expected 43)')
cilia = df[df.name == 'WP_CILIOPATHIES'].iloc[0]
print(f'  WP_CILIOPATHIES degree_norm: {cilia.degree_norm:.3e} (expected 8.509e-06)')
print(f'  WP_CILIOPATHIES null_z: {cilia.null_z:.2f} (expected 41.19)')
print(f'  WP_CILIOPATHIES empirical_q: {cilia.empirical_q:.4f} (expected 0.0058)')
print(f'  WP_CILIOPATHIES member_mean_null_z: {cilia.member_mean_null_z:.3f} (expected 1.385)')
print(f'  WP_CILIOPATHIES member_mean_p: {cilia.member_mean_p:.4f} (expected 0.0569)')
"

python3 -c "
import pandas as pd

df = pd.read_csv('results/kf_nbl/pathway_scores_standard.csv')
df = df.sort_values('degree_norm', ascending=False)
df.index = range(1, len(df)+1)
cilia_rank = df[df.name == 'WP_CILIOPATHIES'].index[0]
print('=== KF-NBL RESULTS ===')
print(f'  WP_CILIOPATHIES rank: {cilia_rank} (expected 3)')
cilia = df[df.name == 'WP_CILIOPATHIES'].iloc[0]
print(f'  WP_CILIOPATHIES degree_norm: {cilia.degree_norm:.3e} (expected 4.241e-06)')
print(f'  WP_CILIOPATHIES null_z: {cilia.null_z:.2f} (expected 18.95)')
print(f'  WP_CILIOPATHIES empirical_q: {cilia.empirical_q:.4f} (expected 0.0116)')
print(f'  WP_CILIOPATHIES member_mean_null_z: {cilia.member_mean_null_z:.3f} (expected 2.433)')
print(f'  WP_CILIOPATHIES member_mean_p: {cilia.member_mean_p:.4f} (expected 0.0030)')

import json
rn = json.load(open('results/kf_nbl/resampling_summary.json'))
pn = rn['primary_run']
print()
print('=== KF-NBL resampling primary run ===')
p10n = pn['bifo_p10']; apn = pn['bifo_ap']
sf10n = pn['sf_p10']; sfapn = pn['sf_ap']
print(f'  P@10={p10n} (expected 0.10)')
print(f'  AP={apn:.3f} (expected 0.044)')
print(f'  Fisher P@10={sf10n} (expected 0.20)')
print(f'  Fisher AP={sfapn:.3f} (expected 0.099)')
"
```

---

## Minimal Test Run

To verify the pipeline end-to-end on a small example (no HPC required):

```bash
cd examples/minimal_test
bash run_test.sh
```

**Expected runtime:** ~5 minutes on a laptop. Produces pathway scores and benchmark metrics matching the curated CHD benchmark at reduced scale.

---

## Frozen Parameter Registry

All analysis parameters are documented in `BENCHMARK_MANIFEST.md`. Key settings:

| Parameter | Value |
|-----------|-------|
| PPR alpha | 0.5 |
| PPR tolerance | 1e-10 |
| PPR max iterations | 500 |
| bifo_mapping.yaml version | v0.7.1 |
| min_members filter | 8 |
| max_members filter | 300 |
| Null model (pathway-node) | membership-rewiring |
| Null model (member-level) | stratified gene set permutation |
| Null permutations | 1000 |
| Null random seed | 42 |
| Member null bins | 10 degree × 10 membership count (log-quantile) |
| Member null sampling | without replacement within null set |

See `BENCHMARK_MANIFEST.md` for complete parameter documentation and expected output checksums.

---

## Generating Summary Output Files

After running `score_pathways.py`, use `summarize_results.py` to produce two output files for each cohort or benchmark run:

```bash
python pipeline/summarize_results.py \
    --scores        results/kf_chd/pathway_scores_standard.csv \
    --seeds         data/cohorts/chd/kf_chd_seeds.txt \
    --reference     data/cohorts/chd/kf_chd_cilia_reference.txt \
    --cohort-name   "KF-CHD" \
    --disease       "congenital heart disease" \
    --n-probands    697 \
    --outdir        results/kf_chd/
```

**Output 1: `pathway_results_summary.tsv`**
Machine-readable table of all scored pathways with clean column names. Suitable for direct loading in R, Python, or Excel. One row per pathway, sorted by degree_norm rank. Includes rank, degree_norm, null_z, empirical_q, null_calibrated, member_mean_null_z, and in_reference columns.

**Output 2: `pathway_results_llm.md`**
Structured markdown document for input to a large language model (LLM). Contains a role instruction, biological context, column interpretation guide, top-50 results table, seed gene list, and suggested questions. Paste into ChatGPT, Claude, or any LLM to discuss the biological meaning of results without prior knowledge of BIFO.

For the KF-CHD and KF-NBL cohorts, the seed files are:
- `data/cohorts/chd/kf_chd_seeds.txt`
- `data/cohorts/nbl/kf_nbl_seeds.txt`

Run `python pipeline/summarize_results.py --help` for all options.
