# BENCHMARK_MANIFEST.md
# BIFO Benchmark v1.0 — Frozen Parameter Registry
# Benchmark run: V5 | Date frozen: April 2026
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
| yaml_version | v0.7.1 | `config/bifo_ddkg_mapping.yaml` |
| gene_universe | all C-prefixed nodes (~13,000) | `pipeline/baseline_enrichment.py` |

## Expected output metrics (benchmark run V5)

### Three-arm pathway ablation (CHD curated, 550 pathways, 18-pathway reference)

| Arm | P@10 | Enrich@10 | Mean rank | Rank improvement |
|-----|------|-----------|-----------|-----------------|
| Full (BIFO conditioned) | 0.700 | 21.4× | 113 | +99.1 |
| Ablation (no bridge edges) | 0.600 | 18.3× | 111 | −11.2 |
| Mechanistic-only | 0.000 | 0.0× | 177 | uninterpretable |

### Baseline comparison (CHD curated)

| Method | P@10 | NDCG@10 | Avg. Prec. | Mean rank |
|--------|------|---------|-----------|-----------|
| Degree overlap (B0) | 0.400 | 0.492 | 0.342 | 84 |
| Seed-only Fisher (B1) | 0.300 | 0.215 | 0.156 | 120 |
| Neighborhood Fisher (B2) | 0.000 | 0.000 | 0.037 | 243 |
| Raw PPR GSEA (B3) | 0.100 | 0.220 | 0.117 | 162 |
| Conditioned PPR GSEA (B3b) | 0.100 | 0.085 | 0.114 | 110 |
| BIFO full-arm (B4) | 0.700 | 0.757 | 0.403 | 113 |

### Four-arm gene-level recovery (CHD curated, 34,523 nodes)

| Arm | AUROC | AUPRC | Entropy | Nonzero (%) |
|-----|-------|-------|---------|-------------|
| Raw | 1.000 | 0.2215 | 5.728 | 100.0 |
| Conditioned (BIFO) | 1.000 | 0.1923 | 5.222 | 32.8 |
| Ablation | 1.000 | 0.2215 | 4.939 | 19.5 |
| Random sparsification | 1.000 | 0.2173 | 5.590 | 64.6 |

### Exhaustive resampling (3,003 CHD splits)

| Metric | Mean | SD | Min | Median | Max |
|--------|------|----|-----|--------|-----|
| BIFO P@10 | 0.445 | 0.123 | 0.100 | 0.500 | 0.700 |
| Rank improvement | +93.0 | 16.8 | +28.4 | +93.5 | +139.1 |
| Splits with positive rank improvement | 3003/3003 (100%) | — | — | — | — |
