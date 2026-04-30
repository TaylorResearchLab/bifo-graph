# BENCHMARK_MANIFEST.md
# BIFO Benchmark — Frozen Parameter Registry
# Date frozen: April 2026
# This file is the reproducibility contract for the manuscript benchmark.
# If any output metric differs from the values below, the run is non-reproducible.

## Frozen parameters

| Parameter | Value | Set in |
|-----------|-------|--------|
| alpha (PPR restart) | 0.5 | `run_conditioning.sh` |
| pathway_min_members | 8 | `score_pathways.py` |
| pathway_max_members | 300 | `score_pathways.py` |
| name_pattern_exclusions | `_Q2, _Q3, _Q4, _Q5, _Q6, MIR` | `score_pathways.py` |
| score_variant | `degree_norm` | `score_pathways.py` |
| c4_random_seed | 42 | `data/benchmark/c4_*_seed_nodes.txt` |
| convergence_tolerance | 1e-10 | `pipeline/bifo_conditioning.py` |
| max_ppr_iterations | 500 | `pipeline/bifo_conditioning.py` |
| yaml_version | v0.7.1 | `config/bifo_mapping_ddkg.yaml` |
| gene_universe | all C-prefixed nodes (~13,000) | `pipeline/baseline_enrichment.py` |
| pathway_direction | unidirectional gene→pathway (BIFO spec v0.02) | `config/bifo_mapping_ddkg.yaml` |

## Expected output metrics

### Conditioning (CHD curated, 174,352 merged input edges)

| Metric | Value |
|--------|-------|
| Total kept edges | 105,192 |
| Propagating edges (full arm) | 57,005 |
| Non-propagating retained | 48,187 |
| Pathway Contribution propagating (gene→PW) | 43,698 (76.7% of propagating) |
| Conditioned AUPRC | 0.1947 |
| Conditioned entropy | 4.934 |

### Three-arm pathway ablation (CHD curated, 550 pathways, 18-pathway reference)

| Arm | Prop. edges | P@10 | Enrich@10 | Mean rank | Rank improvement |
|-----|-------------|------|-----------|-----------|-----------------|
| Full (BIFO conditioned) | 57,005 | 0.700 | 21.4× | 86.6 | +125.4 |
| Ablation (no bridge edges) | 14,413 | 0.600 | 18.3× | 86.3 | +13.5 |
| Mechanistic-only | 9,710 | 0.000 | 0.0× | 177.3 | uninterpretable |

### Baseline comparison (CHD curated)

| Method | P@10 | NDCG@10 | Avg. Prec. | Mean rank |
|--------|------|---------|-----------|-----------|
| Degree overlap (B0) | 0.400 | 0.492 | 0.343 | 83 |
| Seed-only Fisher (B1) | 0.300 | 0.215 | 0.156 | 120 |
| Neighborhood Fisher (B2) | 0.000 | 0.000 | 0.037 | 243 |
| Raw PPR GSEA (B3) | 0.100 | 0.220 | 0.117 | 162 |
| Conditioned PPR GSEA (B3b) | 0.100 | 0.220 | 0.088 | 260 |
| BIFO full-arm (B4) | 0.700 | 0.753 | 0.408 | 87 |

### Four-arm gene-level recovery (CHD curated, 34,523 nodes)

| Arm | AUROC | AUPRC | Entropy |
|-----|-------|-------|---------|
| Raw | 1.000 | 0.2215 | 5.728 |
| Conditioned (BIFO) | 1.000 | 0.1947 | 4.934 |
| Ablation | 1.000 | 0.2215 | 4.919 |
| Random sparsification | 1.000 | 0.2173 | 5.590 |

### Exhaustive resampling (3,003 CHD splits)

| Metric | Mean | SD | Min | Median | Max |
|--------|------|----|-----|--------|-----|
| BIFO P@10 | 0.430 | 0.121 | 0.000 | 0.400 | 0.700 |
| Rank improvement | +118.3 | 16.8 | — | +118.4 | — |
| P@10 ≥ 0.30 | 2835/3003 (94.4%) | — | — | — | — |
| P@10 ≥ 0.50 | 1314/3003 (43.8%) | — | — | — | — |
| Splits with positive rank improvement | 3003/3003 (100%) | — | — | — | — |
| Primary split rank improvement | +128.9 | — | — | — | — |
| Primary split P@10 | 0.700 | — | — | — | — |

## Empirical Null Model (Membership Rewiring)

| Parameter | Value |
|-----------|-------|
| null_type | membership-rewiring |
| n_permutations | 1000 |
| n_swaps_multiplier | 10 |
| ppr_alpha | 0.5 |
| random_seed | 42 |
| min_members | 8 |
| max_members | 300 |

### Benchmark Null Results (CHD curated, N=1000 permutations)

| Metric | Value |
|--------|-------|
| Pathways tested | 550 |
| Significant rewiring null (q<0.05) | 16 |
| Significant rewiring null (q<0.10) | 28 |
| Significant member-level null (q<0.05) | 17 |
| BRUNEAU_SEPTATION_VENTRICULAR null_z | 15.6 |
| BRUNEAU_SEPTATION_VENTRICULAR q | 0.034 |
| WP_HEART_DEVELOPMENT null_z | 16.3 |
| WP_HEART_DEVELOPMENT q | 0.034 |
| WP_CARDIAC_PROGENITOR_DIFFERENTIATION null_z | 12.3 |

