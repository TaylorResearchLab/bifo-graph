# bifo-graph

**Biological Information Flow Ontology (BIFO) — Graph Analysis Pipeline**

Code and data repository for:

> Taylor DM, Mohseni Ahooyi T, Stear B, Zhang Y, Lahiri A, Simmons A, Callahan TJ, Silverstein JC.
> *BIFO: A Biological Information Flow Ontology for Knowledge Graph-Directed Pathway Analysis of Rare Variant Cohort Data.*
> Manuscript in preparation, 2026.

BIFO defines admissible biological information flow on knowledge graphs, transforming them into directed propagation substrates for mechanistically grounded biological inference. BIFO operates on propagated signal rather than set membership, enabling inference over biological communication structure rather than direct overlap alone. Applied to variant-derived gene lists, it recovers ranked pathway hypotheses using personalized PageRank (PPR) propagation over the Data Distillery Knowledge Graph (DDKG).

BIFO complements standard enrichment methods: Fisher enrichment identifies direct gene–pathway overlap, while BIFO identifies propagated signal across biological relationships. When both methods identify the same top pathway independently, their convergence provides strong validation of the biological signal.

---

## Repository structure

```
bifo-graph/
├── pipeline/                    Core Python analysis scripts
│   ├── bifo_conditioning.py     BIFO edge conditioning + PPR propagation (primary)
│   ├── score_pathways.py        Pathway scoring + empirical null models (degree_norm, member_mean, null significance)
│   ├── baseline_enrichment.py   Enrichment baselines (Fisher, GSEA, degree overlap)
│   ├── chd_resampling_exhaustive.py  Exhaustive CHD seed-split resampling (3,003 splits)
│   ├── kf_resampling.py         Bootstrap resampling for KF cohort analyses
│   ├── seed_cui_lookup.py       Map gene symbols / HGNC IDs to UMLS CUIs
│   ├── generate_export_cypher.py  Generate Neo4j export Cypher from seed list
│   ├── build_cilia_reference.py  Build cilia reference pathway set
│   └── clean_cypher_output.py   Clean cypher-shell CSV output artifacts
│
├── scripts/                     Shell wrappers — one per pipeline stage
│   ├── run_full_pipeline.sh     End-to-end pipeline (all stages, one cohort)
│   ├── setup_workspace.sh       Create output directory structure
│   ├── run_kf_chd_export.sh     Stage 1: Neo4j export — KF-CHD cohort
│   ├── run_kf_nbl_export.sh     Stage 1: Neo4j export — KF-NBL cohort
│   ├── clean_files.sh           Stage 2.1: Clean cypher-shell output
│   ├── merge_files.sh           Stage 2.3: Merge edge CSV files
│   ├── run_seed_lookup.sh       Stage 2.4: Map gene symbols to CUIs
│   ├── run_conditioning.sh      Stage 3: BIFO conditioning + PPR propagation
│   ├── run_scoring.sh           Stage 4: Pathway scoring
│   ├── run_baseline.sh          Stage 5: Baseline enrichment comparison
│   ├── run_resampling.sh        Stage 6: Bootstrap resampling
│   └── build_cilia_ref.sh       Build cilia reference set (one-time setup)
│
├── config/
│   └── bifo_mapping.yaml   BIFO edge admissibility rules (v0.7.1)
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
│       │   └── kf_chd_cilia_reference.txt  MSigDB cilia pathway reference set (17 CUIs)
│       │
│       └── nbl/                 KF-NBL variant gene seeds and references
│           ├── kf_nbl_seeds_maf001.txt     Primary seeds (MAF≤0.001, n≥1; 1,406 genes)
│           ├── kf_nbl_seeds_maf001_n2.txt  Carrier-filtered seeds (n≥2; 401 genes)
│           ├── kf_nbl_seeds_maf001_n3.txt  Carrier-filtered seeds (n≥3; 147 genes)
│           ├── kf_nbl_seeds_maf01.txt      MAF≤0.01 seeds (broader filter)
│           ├── kf_nbl_seeds.txt            Original unfiltered seed list
│           ├── kf_nbl_seed_cuis.txt        CUI-resolved seed list (pipeline input)
│           └── kf_nbl_cilia_reference.txt  MSigDB cilia pathway reference set (17 CUIs)
│
├── examples/
│   └── minimal_test/            Self-contained 15-gene CHD end-to-end test
│       ├── README.md            Instructions for the minimal test run
│       ├── neo4j_export.cypher  Export query for minimal test graph
│       ├── seed_nodes.txt       10 CHD seed genes (HGNC CUIs)
│       ├── heldout_nodes.txt    5 held-out CHD genes
│       └── run_test.sh          One-command test run script
│
├── BENCHMARK_MANIFEST.md        Frozen parameter registry and expected output metrics
├── REPRODUCE.md                 Exact commands to reproduce all manuscript results
├── requirements.txt             Python dependencies
│
├── data/
│   ├── benchmark/               (existing — see above)
│   └── cohorts/
│       ├── README.md            Cohort data documentation
│       ├── chd/                 (existing — see above)
│       └── nbl/
│           ├── kf_nbl_seeds_maf001.txt     Primary seeds (MAF≤0.001, n≥1; 1,406 genes)
│           ├── kf_nbl_seeds_maf001_n2.txt  Carrier-filtered seeds (n≥2; 401 genes)
│           ├── kf_nbl_seeds_maf001_n3.txt  Carrier-filtered seeds (n≥3; 147 genes)
│           ├── kf_nbl_seeds_maf01.txt      MAF≤0.01 seeds (broader filter)
│           ├── kf_nbl_seeds.txt            Original unfiltered seed list
│           ├── kf_nbl_seed_cuis.txt        CUI-resolved seed list (pipeline input)
│           └── kf_nbl_cilia_reference.txt  MSigDB cilia pathway reference set (17 CUIs)
│
└── results/                     Pre-computed results (shipped with repo)
    ├── chd_benchmark/
    │   ├── results_full.json               Full-arm PPR results (gene-level metrics)
    │   ├── results_full_scores_cond.npy    Full-arm conditioned PPR scores
    │   ├── results_full_scores_raw.npy     Full-arm raw PPR scores
    │   ├── results_full_kept_edges.csv.zip Full-arm kept edges
    │   ├── results_ablation.json           Ablation-arm PPR results
    │   ├── results_ablation_scores_cond.npy
    │   ├── results_ablation_scores_raw.npy
    │   ├── results_ablation_kept_edges.csv.zip
    │   ├── results_mech.json               Mechanistic-arm PPR results
    │   ├── results_mech_scores_cond.npy
    │   ├── results_mech_scores_raw.npy
    │   ├── results_mech_kept_edges.csv.zip
    │   ├── results_node_index.json         Node CUI → index mapping
    │   ├── pathway_metrics_full.json       Full-arm pathway scoring metrics
    │   ├── pathway_metrics_ablation.json   Ablation-arm pathway scoring metrics
    │   ├── pathway_metrics_mech.json       Mechanistic-arm pathway scoring metrics
    │   ├── baseline_comparison.csv         Baseline enrichment comparison
    │   ├── baseline_comparison.json
    │   ├── resampling_summary.json         Exhaustive 3,003-split resampling summary
    │   ├── resampling_results.csv          Per-split resampling results
    │   ├── test_structural_null.csv        Benchmark null model scores (both nulls)
    │   ├── c4_notch/                       C4/Notch pathway-split control results
    │   │   ├── pathway_metrics.json
    │   │   ├── baseline_comparison.json
    │   │   ├── results.json
    │   │   ├── results_kept_edges.csv.zip
    │   │   ├── results_node_index.json
    │   │   ├── results_scores_cond.npy
    │   │   └── results_scores_raw.npy
    │   └── c4_mapk/                        C4/MAPK pathway-split control results
    │       ├── pathway_metrics.json
    │       ├── baseline_comparison.json
    │       ├── results.json
    │       ├── results_kept_edges.csv.zip
    │       ├── results_node_index.json
    │       ├── results_scores_cond.npy
    │       └── results_scores_raw.npy
    ├── kf_chd/
    │   ├── results.json                    PPR propagation results
    │   ├── results_node_index.json         Node CUI → index mapping
    │   ├── results_kept_edges.csv.zip      BIFO-conditioned kept edges
    │   ├── results_scores_cond.npy         Conditioned PPR score vector
    │   ├── results_scores_raw.npy          Raw PPR score vector
    │   ├── results_scores_meta.npy         Metadata-filtered PPR score vector
    │   ├── results_scores_rand.npy         Random sparsification score vector
    │   ├── nodes_clean_noncc.csv.gz        Node table (pipeline input)
    │   ├── edges_all_noncc.csv.gz          Edge table (pipeline input)
    │   ├── pathway_scores_standard.csv     Pathway scores (full universe)
        │   ├── pathway_scores_null.csv         Pathway scores with empirical null results
    │   ├── pathway_metrics_standard.json   Pathway scoring metrics
        │   ├── baseline_comparison.csv         Baseline enrichment comparison
    │   ├── baseline_comparison.json
    │   ├── resampling_summary.json         Bootstrap resampling summary
    │   └── resampling_results.csv          Per-bootstrap resampling results
    └── kf_nbl/
        ├── results.json                    PPR propagation results
        ├── results_node_index.json         Node CUI → index mapping
        ├── results_kept_edges.csv.zip      BIFO-conditioned kept edges
        ├── results_scores_cond.npy         Conditioned PPR score vector
        ├── results_scores_raw.npy          Raw PPR score vector
        ├── results_scores_meta.npy         Metadata-filtered PPR score vector
        ├── results_scores_rand.npy         Random sparsification score vector
        ├── nodes_clean_noncc.csv.gz        Node table (pipeline input)
        ├── edges_all_noncc.csv.gz          Edge table (pipeline input)
        ├── pathway_scores_standard.csv     Pathway scores (full universe)
                ├── pathway_scores_null.csv         Pathway scores with empirical null results
        ├── pathway_metrics_standard.json   Pathway scoring metrics
                ├── baseline_comparison.csv         Baseline enrichment comparison
        ├── baseline_comparison.json
        ├── resampling_summary.json         Bootstrap resampling summary
        └── resampling_results.csv          Per-bootstrap resampling results
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
# 1. Run BIFO conditioning + PPR propagation
python pipeline/bifo_conditioning.py \
  --nodes data/benchmark/chd_curated_nodes.csv \
  --edges data/benchmark/chd_curated_edges_raw.csv.zip \
  --mapping config/bifo_mapping.yaml \
  --seed-nodes data/benchmark/chd_seed_nodes.txt \
  --heldout-nodes data/benchmark/chd_heldout_nodes.txt \
  --out-json results/chd_benchmark/results.json

# Outputs written alongside results.json:
#   results_scores_cond.npy, results_scores_raw.npy,
#   results_node_index.json, results_kept_edges.csv

# 2. Score pathways
python pipeline/score_pathways.py \
  --nodes data/benchmark/chd_curated_nodes.csv \
  --edges-raw data/benchmark/chd_curated_edges_raw.csv.zip \
  --edges-conditioned results/chd_benchmark/results_kept_edges.csv \
  --scores-cond results/chd_benchmark/results_scores_cond.npy \
  --scores-raw results/chd_benchmark/results_scores_raw.npy \
  --node-index results/chd_benchmark/results_node_index.json \
  --seed-nodes data/benchmark/chd_seed_nodes.txt \
  --chd-pathways data/benchmark/chd_pathway_reference.txt \
  --out-csv results/chd_benchmark/pathway_scores_standard.csv \
  --out-json results/chd_benchmark/pathway_metrics_standard.json

# Note on Quick Start vs paper-faithful run:
# This Quick Start uses edges_raw.csv for simplicity and is a sanity-check
# run only. The paper's primary full-arm analysis uses edges_merged.csv
# (raw mechanistic edges + pathway membership edges combined). See REPRODUCE.md
# for the exact commands used in the manuscript.

# 3. Run baselines
python pipeline/baseline_enrichment.py \
  --edges-merged data/benchmark/chd_curated_edges_raw.csv.zip \
  --node-index results/chd_benchmark/results_node_index.json \
  --scores-raw results/chd_benchmark/results_scores_raw.npy \
  --scores-cond results/chd_benchmark/results_scores_cond.npy \
  --bifo-scores results/chd_benchmark/pathway_scores_standard.csv \
  --chd-pathways data/benchmark/chd_pathway_reference.txt \
  --seed-nodes data/benchmark/chd_seed_nodes.txt \
  --out-csv results/chd_benchmark/baseline_comparison.csv \
  --out-json results/chd_benchmark/baseline_comparison.json

# Expected outputs after steps 1-2:
#   pathway_scores_standard.csv — 550 ranked pathways
#   BRUNEAU_SEPTATION_VENTRICULAR rank 1
#   P@10 = 0.70, mean rank improvement = +125.4

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

## Mapping to manuscript sections

| Methods section | Script / component |
|-----------------|-------------------|
| §1 Knowledge graph source | `cypher/` export queries, `pipeline/clean_cypher_output.py` |
| §2 BIFO conditioning | `pipeline/bifo_conditioning.py`, `config/bifo_mapping.yaml` |
| §3 Personalized PageRank | `pipeline/bifo_conditioning.py` (PPR loop) |
| §4 Pathway scoring | `pipeline/score_pathways.py` |
| §5 Benchmark design | `data/benchmark/`, `BENCHMARK_MANIFEST.md` |
| §6 Baseline enrichment | `pipeline/baseline_enrichment.py` |
| §9 Exhaustive resampling | `pipeline/chd_resampling_exhaustive.py` |
| §10 KF cohort analysis | `scripts/run_full_pipeline.sh`, `pipeline/kf_resampling.py` |
| §8.4 Empirical null models | `pipeline/score_pathways.py` (`--n-permutations`, `--null-type`) |

---

## Pipeline stages

| Stage | Script | Input | Output |
|-------|--------|-------|--------|
| 1 | `run_kf_{cohort}_export.sh` | Neo4j DDKG instance | `edges_raw.csv`, `nodes.csv`, `pathway_membership_edges.csv` |
| 2.1 | `clean_files.sh` | Raw cypher-shell CSVs | Cleaned CSVs |
| 2.3 | `merge_files.sh` | Edge CSVs | `edges_merged.csv` |
| 2.4 | `run_seed_lookup.sh` | Gene symbol list | `seed_cuis.txt` |
| 3 | `run_conditioning.sh` | `edges_merged.csv`, `nodes.csv`, seeds | `results_scores_*.npy`, `results_kept_edges.csv` |
| 4 | `run_scoring.sh` | Score vectors, pathway membership | `pathway_scores_*.csv`, `pathway_metrics_*.json` |
| 5 | `run_baseline.sh` | Edge files, seeds, pathway reference | `baseline_comparison.csv`, `baseline_comparison.json` |
| 6 | `run_resampling.sh` | Kept edges, pathway reference | `resampling_summary.json` |

---

## Empirical null models

`score_pathways.py` implements two complementary empirical null frameworks:

**Pathway-node null (membership rewiring)**
Degree-preserving rewiring of gene→pathway bridge edges with PPR reruns. Tests whether a pathway's concept node receives more propagated mass than expected under randomized membership. Calibration depends on graph composition: valid when non-bridge edges provide sufficient routing constraint (benchmark: 76.7% bridge; KF-CHD: 41.4% bridge; KF-NBL: similar). Under the unidirectional pipeline all three analyses have well-calibrated rewiring nulls. Percentages derived from propagating edge counts in each conditioned graph.

**Member-level null (stratified gene set permutation)**
Matches genes on structural features only (conditioned graph degree and pathway membership count, both log-binned) and tests whether pathway member genes carry disproportionate propagated signal relative to matched random gene sets. Operates on the fixed propagated score vector — no PPR reruns required. Less sensitive to bridge edge fraction than the pathway-node rewiring null, and empirically stable across seed sizes. Tests whether propagated signal concentrates within pathway member genes, not whether signal reaches the pathway node itself.

Both nulls run in a single invocation when `--null-type membership-rewiring` is specified:

```bash
python pipeline/score_pathways.py \
  ... \
  --n-permutations 1000 \
  --null-type membership-rewiring \
  --n-cores 120
```

Output columns: `empirical_q`, `null_z` (pathway-node null); `member_mean_q`, `member_mean_null_z` (member-level null).

The pathway-node null and member-level null test distinct quantities: the former evaluates pathway-node signal accumulation, while the latter evaluates concentration of propagated signal within pathway member genes.

See Methods §8.4 and REPRODUCE.md for full validation results.

---

## Benchmark freeze

Default frozen parameters used in all manuscript analyses:

| Parameter | Value |
|-----------|-------|
| PPR alpha (α) | 0.5 |
| PPR tolerance | 1×10⁻¹⁰ |
| PPR max iterations | 500 |
| min_members filter | 8 |
| max_members filter | 300 |
| Null permutations | 1000 |
| Null random seed | 42 |

Default pathway name exclusions: `*_Q2`–`*_Q6`, `MIR*` (suppresses quartile gene sets and miRNA targets).

Frozen benchmark parameters and expected output metrics are documented in
`BENCHMARK_MANIFEST.md` at the repository root. This file serves as the
reproducibility contract for the manuscript: if any output metric differs
from the manifest values, the run is considered non-reproducible.

Pipeline scripts are in `pipeline/`, the YAML mapping (v0.7.1) is in
`config/`, and seed, heldout, and reference files for all benchmark cohorts
are in `data/benchmark/`.

---

## Key configuration: `config/bifo_mapping.yaml`

The YAML file encodes the BIFO flow class definitions (v0.7.1):
- **251** predicate-to-flow class entries
- **96** explicit non-flow (excluded) predicate designations
- **46** observational edge definitions
- **5** classification tiers: `mechanistic`, `weak_mechanistic_or_observational`,
  `observational`, `contextual_constraint`, `nonpropagating_context` (pathway→gene direction excluded per BIFO spec v0.02)

This file encodes the operational instantiation of BIFO for the DDKG and is
the primary artifact controlling propagation behavior in this analysis.

**Note on two-layer architecture:** The observation that mechanistic edges alone
cannot reach pathway nodes reflects the structure of the DDKG export used here,
in which gene-to-pathway membership edges form a structurally distinct layer.
This property may vary with different graph constructions or vocabulary selections. Modifying
it changes which edges are admissible for propagation and will alter all
downstream results.

---

## Exact reproduction commands

See **REPRODUCE.md** at the repository root for exact CLI commands to reproduce
every result reported in the manuscript, including expected output values for
all three analyses (benchmark, KF-CHD, KF-NBL).

Approximate runtimes (120–192 cores, HPC):
- Benchmark scoring + 1000-permutation null: ~3 minutes
- KF-CHD scoring + 1000-permutation null: ~35 minutes
- KF-NBL scoring + 1000-permutation null: ~55 minutes

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
- KF-CHD and KF-NBL variant gene seed lists (including CUI-resolved seed files)
- Pre-computed KF-CHD and KF-NBL conditioning results, score vectors, pathway scores, and null model results
- Exhaustive CHD resampling summary (3,003 splits) and KF cohort bootstrap resampling summaries
- Benchmark and cohort null model scored outputs (`test_structural_null.csv`, `pathway_scores_null.csv`)


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
