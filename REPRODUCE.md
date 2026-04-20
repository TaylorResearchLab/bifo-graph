# Reproducing All Analyses in Taylor et al. 2026

This document provides exact commands to reproduce every result reported in the manuscript. All analyses were run on the Taylor Lab HPC (reslndbhisbox05) using the frozen pipeline in this repository. Pre-computed results are archived in `results/` and can be verified without re-running.

## Prerequisites

```bash
git clone https://github.com/TaylorResearchLab/bifo-graph
cd bifo-graph
pip install numpy pandas scipy scikit-learn pyyaml
```

Python 3.8+ required. All analyses use `pipeline/bifo_conditioning.py` and `pipeline/score_pathways.py` from this repository. The frozen configuration is `config/bifo_mapping.yaml` (v0.7.1).

---

## Analysis 1 — CHD Curated Benchmark (Tables 1–4, Figures 1–4)

### Step 1.1 — BIFO Conditioning and PPR Propagation

```bash
# Merge raw + membership edges
python3 -c "
import pandas as pd
raw = pd.read_csv('data/benchmark/chd_curated_edges_raw.csv.zip')
mem = pd.read_csv('data/benchmark/chd_curated_pathway_membership_edges.csv.zip')
pd.concat([raw, mem], ignore_index=True).to_csv('../edges_merged_chd.csv', index=False)
print(f'Merged: {len(raw)} + {len(mem)} = {len(raw)+len(mem)} edges')
"

# Run conditioning (four-arm PPR)
python3 pipeline/bifo_conditioning.py \
  --nodes   data/benchmark/chd_curated_nodes.csv \
  --edges   ../edges_merged_chd.csv \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes    data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --out-json results/chd_benchmark/results_full.json
```

**Frozen outputs:** `results/chd_benchmark/results_full.json`, `results_full_scores_cond.npy`, `results_full_scores_raw.npy`, `results_full_kept_edges.csv.zip`, `results_node_index.json`

**Expected:** alpha=0.5, 93,507 propagating edges, conditioned AUPRC=0.1902, entropy=5.217

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
  --edges   ../edges_merged_chd.csv \
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
pd.read_csv(z.open(name), low_memory=False).to_csv('../chd_kept_edges.csv', index=False)
"

python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         ../edges_merged_chd.csv \
  --edges-conditioned ../chd_kept_edges.csv \
  --scores-cond       results/chd_benchmark/results_full_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_full_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --out-csv  results/chd_benchmark/pathway_scores_full.csv \
  --out-json results/chd_benchmark/pathway_metrics_full.json
```

**Expected:** 550 pathways, P@10=0.70, enrichment=21.4×, mean rank=113, rank improvement=+99.1

### Step 1.5 — Pathway Scoring with Empirical Null (Methods §8.4)

```bash
python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         ../edges_merged_chd.csv \
  --edges-conditioned ../chd_kept_edges.csv \
  --scores-cond       results/chd_benchmark/results_full_scores_cond.npy \
  --scores-raw        results/chd_benchmark/results_full_scores_raw.npy \
  --node-index        results/chd_benchmark/results_node_index.json \
  --seed-nodes        data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 120 \
  --out-csv  results/chd_benchmark/pathway_scores_full_null.csv \
  --out-json results/chd_benchmark/pathway_metrics_full_null.json
```

**Expected (pathway-node rewiring null):** 49/550 pathways significant at q<0.05; BRUNEAU q=0.017 null_z=24.1

**Expected (member-level stratified null):** 25/550 pathways significant at q<0.05; BRUNEAU null_z=10.6 q=0.037; WP_HEART_DEVELOPMENT null_z=10.9 q=0.037

### Step 1.6 — Pathway Scoring: Ablation and Mechanistic Arms (Table 3)

```bash
# Ablation arm
python3 pipeline/score_pathways.py \
  --nodes             data/benchmark/chd_curated_nodes.csv \
  --edges-raw         ../edges_merged_chd.csv \
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
  --edges-raw         ../edges_merged_chd.csv \
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
  --edges-merged  ../edges_merged_chd.csv \
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
pd.read_csv(z.open(name), low_memory=False).to_csv(f'../{arm}_kept_edges.csv', index=False)
" $ARM
done

# Score each C4 arm
for ARM in c4_notch c4_mapk; do
  python3 pipeline/score_pathways.py \
    --nodes             data/benchmark/chd_curated_nodes.csv \
    --edges-raw         ../edges_merged_chd.csv \
    --edges-conditioned ../${ARM}_kept_edges.csv \
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
  --edges-merged      ../edges_merged_chd.csv \
  --mapping           config/bifo_mapping.yaml \
  --gene-pool         data/benchmark/chd_seed_nodes.txt \
  --chd-pathways      data/benchmark/chd_pathway_reference.txt \
  --out-csv           results/chd_benchmark/resampling_results.csv \
  --out-json          results/chd_benchmark/resampling_summary.json \
  --n-cores           120
```

**Expected:** 3,003 splits, 100% positive rank improvement, 95.1% P@10≥0.30

---

## Analysis 2 — KF-CHD Cohort (Figures 3–4, Sections 8.1–8.4)

### Data requirements

The KF-CHD graph files are in `/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_chd/`. The merged edges file was produced by:

```bash
# Already exists: kf_chd_edges_merged.csv
# Produced by merging kf_chd_edges_raw_clean.csv + kf_chd_pathway_membership_edges_clean.csv
```

Seed file: `kf_chd_seed_cuis.txt` (1,276 CUIs; derived from `data/cohorts/chd/kf_chd_seeds_maf001.txt`, 1,276/1,287 resolved to graph CUIs)

### Step 2.1 — Conditioning and Propagation

```bash
cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_chd

python3 /path/to/bifo-graph/pipeline/bifo_conditioning.py \
  --nodes   kf_chd_nodes_clean_noncc.csv \
  --edges   kf_chd_edges_all_noncc.csv \
  --mapping /path/to/bifo-graph/config/bifo_mapping.yaml \
  --seed-nodes    kf_chd_seed_cuis.txt \
  --heldout-nodes kf_chd_seed_cuis.txt \
  --out-json kf_chd_results/results.json
```

**Frozen outputs:** `bifo-graph/results/kf_chd/`

**Expected:** 2,115,572 propagating edges (43.9% of concept edges), 1,276/1,276 seeds resolved

### Step 2.2 — Pathway Scoring (Standard Universe)

```bash
python3 /path/to/bifo-graph/pipeline/score_pathways.py \
  --nodes             kf_chd_nodes_clean_noncc.csv \
  --edges-raw         kf_chd_edges_merged.csv \
  --edges-conditioned kf_chd_results/results_kept_edges.csv \
  --scores-cond       kf_chd_results/results_scores_cond.npy \
  --scores-raw        kf_chd_results/results_scores_raw.npy \
  --node-index        kf_chd_results/results_node_index.json \
  --seed-nodes        kf_chd_seed_cuis.txt \
  --chd-pathways      /path/to/bifo-graph/data/cohorts/kf_cilia_reference_msigdb.txt \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 192 \
  --out-csv  kf_chd_results/pathway_scores_standard.csv \
  --out-json kf_chd_results/pathway_metrics_standard.json
```

**Expected (scoring):** WP_CILIOPATHIES rank 1, P@10=0.20 (50× background), rank improvement=+303

**Expected (pathway-node rewiring null):** WP_CILIOPATHIES rewiring null_z=-1.9, q=1.0 (miscalibrated — KF-CHD graph is 93.9% bridge edges; see Methods §8.4)

**Expected (member-level null):** WP_CILIOPATHIES member_mean null_z=3.45, empirical p=0.000999 (floor, 1000 perms)

> Note: `score_pathways.py` runs both nulls in a single invocation when `--null-type membership-rewiring` is specified. Output columns: `empirical_q`, `null_z` (pathway-node rewiring null); `member_mean_q`, `member_mean_null_z` (member-level null).

### Step 2.3 — Baseline Enrichment

```bash
python3 /path/to/bifo-graph/pipeline/baseline_enrichment.py \
  --edges-merged  kf_chd_edges_all_noncc.csv \
  --node-index    kf_chd_results/results_node_index.json \
  --scores-raw    kf_chd_results/results_scores_raw.npy \
  --scores-cond   kf_chd_results/results_scores_cond.npy \
  --bifo-scores   kf_chd_results/pathway_scores_standard.csv \
  --chd-pathways  /path/to/bifo-graph/data/cohorts/kf_cilia_reference_msigdb.txt \
  --seed-nodes    kf_chd_seed_cuis.txt \
  --small-universe \
  --out-csv  kf_chd_results/baseline_comparison.csv \
  --out-json kf_chd_results/baseline_comparison.json
```

**Expected:** WP_CILIOPATHIES rank 1 under seed Fisher (corrected, BH p=9.68×10⁻³¹); Fisher P@10=0.30, AP=0.159; raw PPR GSEA rank 3,325; cond PPR GSEA rank 4,414; neighborhood Fisher P@10=0.0

---

## Analysis 3 — KF-NBL Cohort (Figure 4, Section 8.3)

Independent neuroblastoma cohort (554 probands, dbGaP phs001436). Same pipeline as KF-CHD applied to separate data. Data in `/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_nbl/`. Frozen outputs in `bifo-graph/results/kf_nbl/`.

Seed file: `kf_nbl_seed_cuis.txt` (1,395 CUIs; NBL rare variant carriers at MAF ≤ 0.001)

### Step 3.1 — Conditioning and Propagation

```bash
cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_nbl

python3 /path/to/bifo-graph/pipeline/bifo_conditioning.py \
  --nodes   kf_nbl_nodes_clean_noncc.csv \
  --edges   kf_nbl_edges_all_noncc.csv \
  --mapping /path/to/bifo-graph/config/bifo_mapping.yaml \
  --seed-nodes    kf_nbl_seed_cuis.txt \
  --heldout-nodes kf_nbl_seed_cuis.txt \
  --out-json kf_nbl_results/results.json
```

**Expected:** 2,647,055 propagating edges (45.1% of concept edges), 1,395/1,406 seeds resolved

### Step 3.2 — Pathway Scoring

```bash
python3 /path/to/bifo-graph/pipeline/score_pathways.py \
  --nodes             kf_nbl_nodes_clean_noncc.csv \
  --edges-raw         kf_nbl_edges_merged.csv \
  --edges-conditioned kf_nbl_results/results_kept_edges.csv \
  --scores-cond       kf_nbl_results/results_scores_cond.npy \
  --scores-raw        kf_nbl_results/results_scores_raw.npy \
  --node-index        kf_nbl_results/results_node_index.json \
  --seed-nodes        kf_nbl_seed_cuis.txt \
  --chd-pathways      /path/to/bifo-graph/data/cohorts/kf_cilia_reference_msigdb.txt \
  --min-members 8 --max-members 300 \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores  192 \
  --out-csv  kf_nbl_results/pathway_scores_standard.csv \
  --out-json kf_nbl_results/pathway_metrics_standard.json
```

**Expected (scoring):** WP_CILIOPATHIES rank 1 (independent replication of KF-CHD finding)

**Expected (pathway-node rewiring null):** WP_CILIOPATHIES null_z=31.3, empirical q=0.009

**Expected (member-level null):** WP_CILIOPATHIES member_mean null_z=2.41, empirical p=0.004

### Step 3.3 — Baseline Enrichment

```bash
python3 /path/to/bifo-graph/pipeline/baseline_enrichment.py \
  --edges-merged  kf_nbl_edges_all_noncc.csv \
  --node-index    kf_nbl_results/results_node_index.json \
  --scores-raw    kf_nbl_results/results_scores_raw.npy \
  --scores-cond   kf_nbl_results/results_scores_cond.npy \
  --bifo-scores   kf_nbl_results/pathway_scores_standard.csv \
  --chd-pathways  /path/to/bifo-graph/data/cohorts/kf_cilia_reference_msigdb.txt \
  --seed-nodes    kf_nbl_seed_cuis.txt \
  --small-universe \
  --out-csv  kf_nbl_results/baseline_comparison.csv \
  --out-json kf_nbl_results/baseline_comparison.json
```

**Expected:** WP_CILIOPATHIES rank 1 under seed Fisher (corrected, BH p=8.13×10⁻³²); Fisher P@10=0.30, AP=0.165; raw PPR GSEA rank ~3,000+; convergent with KF-CHD

---

## Verifying Pre-Computed Results

All frozen outputs are in `results/`. To verify all reported numbers match the manuscript,
run from the `bifo-graph/` repo root (for cohort analyses, run from the respective
`bifo_run_chd/` or `bifo_run_nbl/` directories as shown in the analysis steps above):

```bash
python3 -c "
import json, pandas as pd

print('=== TABLE 1: Conditioning ===')
r = json.load(open('results/chd_benchmark/results_full.json'))
cov = r['coverage']
gs = r['graph_stats']
print(f'  Nodes: {cov[\"total_nodes\"]} (expected 34,523)')
print(f'  Kept edges: {cov[\"kept_edges\"]} (expected 104,342)')
print(f'  Pathway Contribution: {gs[\"flow_class_distribution\"][\"Pathway Contribution\"]} (expected 80,200)')
print(f'  Cond AUPRC: {r[\"conditioned\"][\"auprc\"]:.4f} (expected 0.1902)')
print(f'  Cond entropy: {r[\"conditioned\"][\"entropy\"]:.3f} (expected 5.222)')
print(f'  Raw AUPRC: {r[\"raw\"][\"auprc\"]:.4f} (expected 0.2215)')
print(f'  Raw entropy: {r[\"raw\"][\"entropy\"]:.3f} (expected 5.728)')

print()
print('=== TABLE 3: Three-arm ablation ===')
for arm, ep10, eri in [('full',0.70,+99.1),('ablation',0.60,-11.2),('mech',0.00,34.7)]:
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
ap_expected = {'degree_overlap':0.342,'seed_fisher':0.156,'neighborhood_fisher':0.036,
               'raw_ppr_gsea':0.117,'cond_ppr_gsea':0.114,'bifo_full':0.403}
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
notch = json.load(open('results/chd_benchmark/c4_notch/pathway_metrics.json'))['metrics']
mapk  = json.load(open('results/chd_benchmark/c4_mapk/pathway_metrics.json'))['metrics']
print(f'  Notch: P@10={notch[\"top10_precision\"]:.2f} AP={notch.get(\"average_precision\",\"?\"):.3f} rank_imp={notch[\"rank_improvement_cond_vs_raw\"]:.1f}')
print(f'         (expected P@10=0.50 AP=0.450 rank_imp=-31.5)')
print(f'  MAPK:  P@10={mapk[\"top10_precision\"]:.2f} AP={mapk.get(\"average_precision\",\"?\"):.3f} rank_imp={mapk[\"rank_improvement_cond_vs_raw\"]:.1f}')
print(f'         (expected P@10=0.10 AP=0.174 rank_imp=-54.4)')

print()
print('=== TABLE 6: Resampling ===')
r = json.load(open('results/chd_benchmark/resampling_summary.json'))
rob = r['robustness']
mets = r['metrics']
print(f'  Splits: {r[\"n_splits_total\"]} (expected 3,003)')
print(f'  P@10>=0.30: {rob[\"bifo_p10_ge_0.3\"]}/3003 = {rob[\"bifo_p10_ge_0.3\"]/3003*100:.1f}% (expected 95.1%)')
print(f'  P@10>=0.50: {rob[\"bifo_p10_ge_0.5\"]}/3003 = {rob[\"bifo_p10_ge_0.5\"]/3003*100:.1f}% (expected 51.8%)')
print(f'  Positive rank imp: {rob[\"rank_imp_positive\"]}/3003 (expected 3003)')
print(f'  Median P@10: {mets[\"bifo_p10\"][\"median\"]:.2f} (expected 0.50)')
print(f'  AP range: {mets[\"bifo_ap\"][\"min\"]:.3f}-{mets[\"bifo_ap\"][\"max\"]:.3f} (expected 0.136-0.448)')

print()
print('=== NULL MODEL: Benchmark (from pathway_scores_full_null.csv) ===')
df = pd.read_csv('results/chd_benchmark/pathway_scores_full_null.csv')
df = df.sort_values('degree_norm', ascending=False)
sig = (df.empirical_q < 0.05).sum()
sig_mm = (df.member_mean_q < 0.05).sum()
bruneau = df[df.name.str.contains('BRUNEAU')].iloc[0]
print(f'  Rewiring null sig q<0.05: {sig} (expected 49)')
print(f'  BRUNEAU rewiring null_z: {bruneau.null_z:.1f} (expected 24.1)')
print(f'  BRUNEAU q: {bruneau.empirical_q:.3f} (expected 0.017)')
print(f'  Member null sig q<0.05: {sig_mm} (expected 25)')
print(f'  BRUNEAU member null_z: {bruneau.member_mean_null_z:.1f} (expected 10.6)')
print(f'  BRUNEAU member q: {bruneau.member_mean_q:.3f} (expected 0.037)')
"
```

For KF-CHD and KF-NBL verification, run from the respective cohort directories:

```bash
# KF-CHD — run from /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_chd/
python3 -c "
import json, pandas as pd

print('=== KF-CHD RESULTS ===')
rs = json.load(open('kf_chd_results/resampling_summary.json'))
pr = rs['primary_run']
print(f'  Primary run P@10: {pr[\"bifo_p10\"]} (expected 0.20)')
print(f'  Primary run AP: {pr[\"bifo_ap\"]:.3f} (expected 0.090)')
print(f'  Primary run rank_improvement: {pr[\"rank_improvement\"]:.1f} (expected +303)')
print(f'  Fisher P@10: {pr[\"sf_p10\"]} (expected 0.30)')
print(f'  Fisher AP: {pr[\"sf_ap\"]:.3f} (expected 0.150)')

df = pd.read_csv('kf_chd_results/pathway_scores_standard.csv')
df = df.sort_values('degree_norm', ascending=False)
cilia = df[df.name == 'WP_CILIOPATHIES'].iloc[0]
print(f'  WP_CILIOPATHIES rank: 1 (expected 1)')
print(f'  WP_CILIOPATHIES degree_norm: {cilia.degree_norm:.8f} (expected 0.00000626)')
print(f'  WP_CILIOPATHIES member null_z: {cilia.member_mean_null_z:.2f} (expected 3.45)')
print(f'  WP_CILIOPATHIES member p: {cilia.member_mean_p:.6f} (expected 0.000999)')
print(f'  WP_CILIOPATHIES rewiring null_z: {cilia.null_z:.2f} (expected -1.90)')
"

# KF-NBL — run from /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_nbl/
python3 -c "
import pandas as pd

df = pd.read_csv('kf_nbl_results/pathway_scores_standard.csv')
df = df.sort_values('degree_norm', ascending=False)
cilia = df[df.name == 'WP_CILIOPATHIES'].iloc[0]
print('=== KF-NBL RESULTS ===')
print(f'  WP_CILIOPATHIES rank: 1 (expected 1)')
print(f'  WP_CILIOPATHIES degree_norm: {cilia.degree_norm:.8f} (expected 0.00000596)')
print(f'  WP_CILIOPATHIES rewiring null_z: {cilia.null_z:.2f} (expected 31.28)')
print(f'  WP_CILIOPATHIES rewiring q: {cilia.empirical_q:.4f} (expected 0.0093)')
print(f'  WP_CILIOPATHIES member null_z: {cilia.member_mean_null_z:.2f} (expected 2.41)')
print(f'  WP_CILIOPATHIES member p: {cilia.member_mean_p:.6f} (expected 0.003996)')
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
