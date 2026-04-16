# bifo-graph

**Biological Information Flow Ontology (BIFO) — Graph Analysis Pipeline**

Transforms variant-derived gene lists into ranked pathway hypotheses using
personalized PageRank propagation over a DDKG/UBKG-based knowledge graph.

## Repository Structure

```
bifo-graph/
├── pipeline/          Core Python scripts (bifo_conditioning, score_pathways, etc.)
├── scripts/           Shell scripts — one per pipeline stage
├── config/            bifo_ddkg_mapping.yaml (BIFO edge admissibility rules)
├── cypher/            Neo4j export queries (one file per cohort)
├── data/
│   ├── benchmark/     Curated CHD benchmark seed/reference sets (C4 controls)
│   ├── cohorts/
│   │   ├── chd/       KF-CHD variant gene seeds and NCC reference
│   │   └── nbl/       KF-NBL variant gene seeds and NCC reference
│   └── ncc_cilia_pathways/   20 NCC/cilia gene set files (shipped as package data)
├── examples/
│   └── minimal_test/  15-gene CHD test run with README
└── results/
    ├── chd_benchmark/ Exhaustive resampling summary (CHD curated benchmark)
    └── kf_chd/        KF-CHD conditioning results and score vectors
```

## Quick Start

```bash
# Full pipeline for CHD cohort
bash scripts/run_full_pipeline.sh chd neo4j PASSWORD bolt://localhost:7687

# Full pipeline for NBL cohort
bash scripts/run_full_pipeline.sh nbl neo4j PASSWORD bolt://localhost:7687
```

## Requirements

- Python 3.7+
- pandas >= 0.25, numpy >= 1.16, scipy >= 1.2, pyyaml >= 5.0
- cypher-shell (for Neo4j graph export)
- Neo4j 5.x with UBKG/Petagraph database

## Pipeline Stages

| Stage | Script | Purpose |
|-------|--------|---------|
| 1 | `run_kf_{cohort}_export.sh` | Export graph data from Neo4j |
| 2.1 | `clean_files.sh` | Clean cypher-shell CSV output |
| 2.2 | `build_ncc_edges.sh` | Build NCC/cilia membership edges |
| 2.3 | `merge_files.sh` | Merge CSV files (pandas, not cat) |
| 2.4 | `run_seed_lookup.sh` | Map gene symbols to UMLS CUIs |
| 3 | `run_conditioning.sh` | BIFO conditioning + PPR propagation |
| 4 | `run_scoring.sh` | Pathway scoring (standard + NCC) |
| 5 | `run_baseline.sh` | Baseline enrichment comparison |
| 6 | `run_resampling.sh` | Bootstrap stability resampling |

## Citation

*Manuscript in preparation.* Taylor Research Lab, CHOP/UPenn.
