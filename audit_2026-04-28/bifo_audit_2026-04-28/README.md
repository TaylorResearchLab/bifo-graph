# BIFO-PPR Pipeline Audit (2026-04-28)

Audit directory at `/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph/audit_2026-04-28/`.

This directory contains the diagnostic findings, pipeline fixes, and rerun scripts for the BIFO-PPR pipeline audit conducted April 28, 2026. The audit was triggered by a discrepancy between the manuscript-canonical KF-CHD WP_CILIOPATHIES result (rank 43, null_z 48.7) and a fresh student-mode REPRODUCE.md re-run (rank 3, null_z ~18).

## Audit conclusion in one paragraph

The CHD pipeline is broken end-to-end as committed because (1) the cypher generator emits 4 of 5 required query files, missing Query 6 (pathway_member_nodes), and (2) the committed CHD edges file was generated against an obsolete 55-gene curated seed list rather than the 1,276-gene production seed list. NBL reproduces because NBL never had a seed-correction event. Both findings require regenerating CHD's pipeline outputs from Cypher-export forward, against the corrected production seed list, with appropriate provenance pinning. While we are at it, four other cohorts (NBL, t-ALL, CBTN, RSBD) are being run in the same session to produce the dataset needed for the proposed Paper 3 and Paper 4 of the U24 publication track.

## Directory contents

```
audit_2026-04-28/
├── README.md                                        <- this file
├── docs/
│   └── audit_findings.md                           <- detailed findings
├── pipeline_fixes/
│   └── generate_export_cypher.py                   <- patched: adds Query 6
├── new_export_scripts/
│   ├── pre_rerun_setup.sh                          <- Phase 0: env check + snapshots
│   ├── run_cohort_export.sh                        <- generic per-cohort export
│   ├── run_kf_chd_export.sh                        <- thin wrapper for CHD
│   ├── run_kf_nbl_export.sh                        <- thin wrapper for NBL
│   ├── run_kf_tall_export.sh                       <- thin wrapper for t-ALL
│   ├── run_kf_cbtn_export.sh                       <- thin wrapper for CBTN
│   ├── run_kf_rsbd_export.sh                       <- thin wrapper for RSBD
│   ├── run_one_cohort.sh                           <- full pipeline for one cohort
│   └── run_all_cohorts.sh                          <- driver for all five cohorts
└── snapshots/                                      <- populated by pre_rerun_setup.sh
    ├── pre_rerun_md5s_*.txt
    └── git_status_*.txt
```

The actual rerun outputs land in `/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph/_rerun_2026-04-28/{chd,nbl,tall,cbtn,rsbd}/` (separate from the audit directory; these are working artifacts, not source-controlled reference material).

## Quick start (HPC operator)

```bash
# 1. Activate conda environment with Python 3.11
conda activate bifo-spectral

# 2. Navigate to repo
cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph

# 3. Apply the patched cypher generator
cp audit_2026-04-28/pipeline_fixes/generate_export_cypher.py pipeline/generate_export_cypher.py

# 4. Phase 0: pre-rerun setup (env checks, md5 snapshots)
bash audit_2026-04-28/new_export_scripts/pre_rerun_setup.sh

# 5. Phase 1: rerun KF-CHD and KF-NBL (~4-5 hours per cohort, can run in parallel)
bash audit_2026-04-28/new_export_scripts/run_all_cohorts.sh --cohorts "chd nbl"

# (Verify CHD produces rank 3 / null_z ~18, NBL reproduces canonical within ~3%.)

# 6. Phase 2: rerun t-ALL, CBTN, RSBD (~12-15 hours, can parallelize across HPC nodes)
bash audit_2026-04-28/new_export_scripts/run_all_cohorts.sh --cohorts "tall cbtn rsbd"
```

## Per-cohort wall time estimates

| Stage | Per-cohort time |
|---|---|
| Cypher export (5 queries) | ~30-90 min (Q5 dominates; ~1M+ rows) |
| Clean files | ~1 min |
| Merge files | ~1 min |
| Seed CUI lookup | ~1 min |
| BIFO conditioning | ~20-30 min |
| Pathway scoring (1000 perms, all cores) | ~75-120 min |
| **Subtotal: export through scoring** | **~3-5 hours per cohort** |
| Baseline enrichment | ~10-30 min |
| Resampling | ~30-60 min |
| **Total with downstream stages** | **~4-6 hours per cohort** |

If running on multiple HPC nodes in parallel, total wall time for all five cohorts: ~6 hours.

## Running cohorts in parallel across HPC nodes

For parallel execution across multiple HPC nodes (recommended given access), launch one cohort per node:

```bash
# On node 1
ssh node1
cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph
conda activate bifo-spectral
bash audit_2026-04-28/new_export_scripts/run_one_cohort.sh chd

# On node 2 (in parallel)
ssh node2
cd /mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph
conda activate bifo-spectral
bash audit_2026-04-28/new_export_scripts/run_one_cohort.sh nbl

# ...etc for tall, cbtn, rsbd
```

All five cohorts share the same Neo4j instance (single bolt:// address) so parallelization is bounded by Neo4j throughput, not pipeline structure. If Neo4j becomes a bottleneck, serialize the Stage 1 (export) phase but parallelize Stages 2-4.

## Phase 1 verification expectations

After Phase 1 completes, the WP_CILIOPATHIES summary table should show:

| Cohort | Expected rank | Expected null_z | Expected q | Note |
|---|---|---|---|---|
| chd | ~3 | ~18 | ~0.014 | Matches today's REPRODUCE.md run, NOT the canonical (rank 43) |
| nbl | 3 | 18.95 ± 0.5 | 0.012 ± 0.003 | Reproduces canonical within numerical noise |

If CHD comes back at rank 43 / null_z 48 (the canonical), the bug-hypothesis is wrong and we need a different diagnostic. If NBL diverges, that also indicates something else is going on. Either anomaly should pause Phase 2 until investigated.

## Phase 2 verification expectations

After Phase 2 completes, WP_CILIOPATHIES across all five cohorts should be at relatively similar rank/null_z. The curation-bias hypothesis predicts that this is the expected pattern for any P/LP-filtered seed set scored on the same DDKG. The exact magnitudes depend on cohort size (larger cohorts have more seeds, often producing slightly different null structures), but the rank should be in the top-10 across all five.

If WP_CILIOPATHIES ranks dramatically differently between cohorts (e.g., top-3 in CHD/NBL but rank >100 in t-ALL/CBTN/RSBD), that would indicate cohort-specific signal beyond what the curation-bias hypothesis predicts and would need investigation.

## What changes after the rerun

After Phase 1 + 2 complete and outputs are validated:

1. **Manuscript Paper 1**: Update KF-CHD pathway result to new authoritative numbers (likely rank ~3, null_z ~18). Update any figures/tables that referenced canonical numbers. Add a note in REPRODUCE.md about the corrected pipeline.

2. **U24 memo**: Update Section 1 (where-we-are) with new authoritative numbers. The cohort-similarity discovery (Section 2) becomes stronger because the new numbers show CHD and NBL at very similar pathway rankings, consistent with the curation-bias hypothesis.

3. **Paper 3 (proposed)**: Phase 2 produces the five-cohort dataset that powers this paper's empirical analysis. Cross-cohort pathway rank correlations become a central figure.

4. **Paper 4 (proposed)**: Awaits Paper 3's bias-correction methodology. Can begin outline once Phase 2 completes.

## Contact

Run owner: Deanne Taylor (taylordm@chop.edu)
Audit conducted: April 28, 2026
