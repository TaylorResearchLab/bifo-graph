# bifo-graph

**Biological Information Flow Ontology (BIFO) — Graph Analysis Pipeline**

Code and data repository for:

> Taylor DM, Mohseni Ahooyi T, Stear B, Zhang Y, Lahiri A, Simmons A, Callahan TJ, Silverstein JC.
> *BIFO: A Biological Information Flow Ontology for Knowledge Graph-Directed Pathway Analysis of Rare Variant Cohort Data.*
> Manuscript in preparation, 2026.

BIFO transforms heterogeneous biomedical knowledge graphs into constrained propagation substrates by classifying edges into biologically admissible flow classes. Applied to variant-derived gene lists, it recovers ranked pathway hypotheses using personalized PageRank (PPR) propagation over the Data Distillery Knowledge Graph (DDKG).

---

## Repository structure

```
bifo-graph/
├── pipeline/                    Core Python analysis scripts
│   ├── bifo_conditioning.py     BIFO edge conditioning + PPR propagation (primary)
│   ├── score_pathways.py        Pathway scoring (degree_norm and variants)
│   ├── baseline_enrichment.py   Enrichment baselines (Fisher, GSEA, degree overlap)
│   ├── chd_resampling_exhaustive.py  Exhaustive CHD seed-split resampling (3,003 splits)
│   ├── kf_resampling.py         Bootstrap resampling for KF cohort analyses
│   ├── seed_cui_lookup.py       Map gene symbols / HGNC IDs to UMLS CUIs
│   ├── generate_export_cypher.py  Generate Neo4j export Cypher from seed list
│   ├── build_ncc_membership_edges.py  Build NCC/cilia pathway membership edges
│   ├── build_cilia_reference.py  Build cilia reference pathway set
│   └── clean_cypher_output.py   Clean cypher-shell CSV output artifacts
│
├── scripts/                     Shell wrappers — one per pipeline stage
│   ├── run_full_pipeline.sh     End-to-end pipeline (all stages, one cohort)
│   ├── setup_workspace.sh       Create output directory structure
│   ├── run_kf_chd_export.sh     Stage 1: Neo4j export — KF-CHD cohort
│   ├── run_kf_nbl_export.sh     Stage 1: Neo4j export — KF-NBL cohort
│   ├── clean_files.sh           Stage 2.1: Clean cypher-shell output
│   ├── build_ncc_edges.sh       Stage 2.2: Build NCC/cilia membership edges
│   ├── merge_files.sh           Stage 2.3: Merge edge CSV files
│   ├── run_seed_lookup.sh       Stage 2.4: Map gene symbols to CUIs
│   ├── run_conditioning.sh      Stage 3: BIFO conditioning + PPR propagation
│   ├── run_scoring.sh           Stage 4: Pathway scoring
│   ├── run_baseline.sh          Stage 5: Baseline enrichment comparison
│   ├── run_resampling.sh        Stage 6: Bootstrap resampling
│   └── build_cilia_ref.sh       Build cilia reference set (one-time setup)
│
├── config/
│   └── bifo_ddkg_mapping.yaml   BIFO edge admissibility rules (v0.7.1)
│                                251 predicate-to-flow entries, 5 classification tiers
│
├── cypher/                      Neo4j export queries (one file per benchmark/cohort)
│   ├── chd_curated_export_queries.cypher   CHD curated benchmark
│   ├── c4_notch_export_queries.cypher      C4/Notch pathway-split control
│   ├── c4_mapk_export_queries.cypher       C4/MAPK pathway-split control
│   ├── kf_chd_export_queries.cypher        KF-CHD cohort (Kids First)
│   └── kf_nbl_export_queries.cypher        KF-NBL cohort (Kids First)
│
├── data/
│   ├── benchmark/               Curated CHD benchmark and C4 control data
│   │   ├── chd_curated_edges_raw.csv.zip          1-hop CHD neighborhood edges
│   │   ├── chd_curated_pathway_membership_edges.csv.zip  Pathway membership edges
│   │   ├── chd_pathway_reference.txt              18-pathway CHD reference set
│   │   ├── c4_notch_seed_nodes.txt                C4/Notch seeds (30 genes)
│   │   ├── c4_notch_heldout_nodes.txt             C4/Notch heldout (14 genes)
│   │   ├── c4_notch_pathway_reference.txt         C4/Notch reference pathways
│   │   ├── c4_mapk_seed_nodes.txt                 C4/MAPK seeds (63 genes)
│   │   ├── c4_mapk_heldout_nodes.txt              C4/MAPK heldout (28 genes)
│   │   └── c4_mapk_pathway_reference.txt          C4/MAPK reference pathways
│   │
│   └── cohorts/
│       ├── chd/                 KF-CHD variant gene seeds and references
│       │   ├── kf_chd_seeds_maf001.txt     Primary seeds (MAF≤0.001, n≥1; 1,287 genes)
│       │   ├── kf_chd_seeds_maf001_n2.txt  Carrier-filtered seeds (n≥2; 387 genes)
│       │   ├── kf_chd_seeds_maf001_n3.txt  Carrier-filtered seeds (n≥3; 146 genes)
│       │   ├── kf_chd_seeds_maf01.txt      MAF≤0.01 seeds (broader filter)
│       │   ├── kf_chd_seeds.txt            Original unfiltered seed list
│       │   ├── kf_chd_seed_cuis.txt        CUI-resolved seed list (pipeline input)
│       │   ├── kf_chd_seed_nodes.csv.zip   Full seed node metadata
│       │   └── kf_chd_ncc_reference.txt    NCC cilia pathway reference set
│       │
│       └── nbl/                 KF-NBL variant gene seeds and references
│           ├── kf_nbl_seeds_maf001.txt     Primary seeds (MAF≤0.001, n≥1; 1,406 genes)
│           ├── kf_nbl_seeds_maf001_n2.txt  Carrier-filtered seeds (n≥2; 401 genes)
│           ├── kf_nbl_seeds_maf001_n3.txt  Carrier-filtered seeds (n≥3; 147 genes)
│           ├── kf_nbl_seeds_maf01.txt      MAF≤0.01 seeds (broader filter)
│           ├── kf_nbl_seeds.txt            Original unfiltered seed list
│           └── kf_nbl_ncc_reference.txt    NCC cilia pathway reference set
│
├── examples/
│   └── minimal_test/            Self-contained 15-gene CHD end-to-end test
│       ├── README.md            Instructions for the minimal test run
│       ├── neo4j_export.cypher  Export query for minimal test graph
│       ├── seed_nodes.txt       10 CHD seed genes (HGNC CUIs)
│       ├── heldout_nodes.txt    5 held-out CHD genes
│       └── run_test.sh          One-command test run script
│
└── results/                     Pre-computed results (shipped with repo)
    ├── chd_benchmark/
    │   └── chd_resampling_summary.json   Exhaustive 3,003-split resampling summary
    └── kf_chd/
        ├── results.json                  Four-arm PPR metrics (gene-level recovery)
        ├── results_node_index.json       Node CUI → index mapping
        ├── results_kept_edges.csv.zip    BIFO-conditioned kept edges
        ├── results_scores_cond.npy       Conditioned PPR score vector
        ├── results_scores_raw.npy        Raw PPR score vector
        ├── results_scores_meta.npy       Metadata-filtered PPR score vector
        └── results_scores_rand.npy       Random sparsification score vector
```

---

## Quick start

### Prerequisites

```bash
# Python dependencies
pip install pandas>=0.25 numpy>=1.16 scipy>=1.2 pyyaml>=5.0

# Neo4j access (for graph export only)
# cypher-shell must be on PATH, connected to a DDKG/UBKG/Petagraph instance
```

### Reproduce the curated CHD benchmark

The curated benchmark graph data (edges and pathway memberships) is shipped
with the repository in `data/benchmark/`. No Neo4j connection is required to
reproduce the primary benchmark results.

```bash
# 1. Run BIFO conditioning + PPR
python pipeline/bifo_conditioning.py \
  --edges-merged data/benchmark/chd_curated_edges_raw.csv.zip \
  --pathway-edges data/benchmark/chd_curated_pathway_membership_edges.csv.zip \
  --nodes data/benchmark/chd_curated_nodes.csv \
  --seed-nodes data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --config config/bifo_ddkg_mapping.yaml \
  --out-stem results/chd_benchmark/results

# 2. Score pathways
python pipeline/score_pathways.py \
  --scores results/chd_benchmark/results_scores_cond.npy \
  --kept-edges results/chd_benchmark/results_kept_edges.csv \
  --pathway-reference data/benchmark/chd_pathway_reference.txt \
  --out results/chd_benchmark/pathway_scores.csv

# 3. Run baselines
python pipeline/baseline_enrichment.py \
  --edges-merged data/benchmark/chd_curated_edges_raw.csv.zip \
  --kept-edges results/chd_benchmark/results_kept_edges.csv \
  --seed-nodes data/benchmark/chd_seed_nodes.txt \
  --pathway-reference data/benchmark/chd_pathway_reference.txt \
  --out results/chd_benchmark/baseline_comparison.json

# 4. Run exhaustive resampling (3,003 splits)
python pipeline/chd_resampling_exhaustive.py \
  --kept-edges results/chd_benchmark/results_kept_edges.csv \
  --edges-merged data/benchmark/chd_curated_edges_raw.csv.zip \
  --pathway-reference data/benchmark/chd_pathway_reference.txt \
  --out results/chd_benchmark/chd_resampling_results.csv \
  --n-cores 4
```

### Run Kids First cohort analysis (requires Neo4j)

```bash
# Full KF-CHD pipeline (all stages)
# Full KF-CHD pipeline (all stages)
bash scripts/run_full_pipeline.sh chd maf001 0 neo4j PASSWORD <your_bolt_instance_address> 

# Full KF-NBL pipeline
bash scripts/run_full_pipeline.sh nbl maf001 0 neo4j PASSWORD <your_bolt_instance_address>
```

### Minimal end-to-end test (no Neo4j required for scoring)

```bash
cd examples/minimal_test/
bash run_test.sh
```

---

## Pipeline stages

| Stage | Script | Input | Output |
|-------|--------|-------|--------|
| 1 | `run_kf_{cohort}_export.sh` | Neo4j DDKG instance | `edges_raw.csv`, `nodes.csv`, `pathway_membership_edges.csv` |
| 2.1 | `clean_files.sh` | Raw cypher-shell CSVs | Cleaned CSVs |
| 2.2 | `build_ncc_edges.sh` | NCC gene sets | `ncc_membership_edges.csv` |
| 2.3 | `merge_files.sh` | Edge CSVs | `edges_merged.csv` |
| 2.4 | `run_seed_lookup.sh` | Gene symbol list | `seed_cuis.txt` |
| 3 | `run_conditioning.sh` | `edges_merged.csv`, `nodes.csv`, seeds | `results_scores_*.npy`, `results_kept_edges.csv` |
| 4 | `run_scoring.sh` | Score vectors, pathway membership | `pathway_scores_*.csv`, `pathway_metrics_*.json` |
| 5 | `run_baseline.sh` | Edge files, seeds, pathway reference | `baseline_comparison.json` |
| 6 | `run_resampling.sh` | Kept edges, pathway reference | `resampling_summary.json` |

---

## Key configuration: `config/bifo_ddkg_mapping.yaml`

The YAML file encodes the BIFO flow class definitions (v0.7.1):
- **251** predicate-to-flow class entries
- **96** explicit non-flow (excluded) predicate designations
- **46** observational edge definitions
- **5** classification tiers: `mechanistic`, `weak_mechanistic_or_observational`,
  `observational`, `contextual_constraint`, `nonpropagating_context`

This file is the primary scientific artifact of the BIFO framework. Modifying
it changes which edges are admissible for propagation and will alter all
downstream results.

---

## Reproducing paper figures

All figures in the manuscript are generated from outputs of the pipeline scripts
above. The pre-computed KF-CHD conditioning results in `results/kf_chd/` can be
used to reproduce figures without re-running the full export pipeline.

Figure source scripts are available in the companion figure repository
(linked from the paper).

---

## Data availability

**Shipped with this repository:**
- Curated CHD benchmark graph (edges and pathway memberships)
- C4 pathway-split control seeds and references
- KF-CHD and KF-NBL variant gene seed lists
- Pre-computed KF-CHD conditioning results and score vectors
- Exhaustive CHD resampling summary (3,003 splits)

**Requires DDKG access (not shipped):**
- Full KF-CHD and KF-NBL graph exports (815K–880K nodes, 5–6M edges)
- Neo4j DDKG/UBKG/Petagraph instance for re-export

KF cohort variant data: Kids First Data Resource Portal
(phs001138 KF-CHD; phs001436 KF-NBL). Access requires dbGaP authorization.

---

## Citation

```bibtex
@article{taylor2026bifo,
  title   = {{BIFO}: A Biological Information Flow Ontology for
             Knowledge Graph-Directed Pathway Analysis of Rare Variant Cohort Data},
  author  = {Taylor, Deanne M. and Mohseni Ahooyi, Taha and Stear, Benjamin
             and Zhang, Yuanchao and Lahiri, Aditya and Simmons, Alan
             and Callahan, Tiffany J. and Silverstein, Jonathan C.},
  journal = {Manuscript in preparation},
  year    = {2026}
}
```

---

## License

Code: MIT License. See `LICENSE`.
Data (benchmark gene sets, pathway references): CC BY 4.0.
