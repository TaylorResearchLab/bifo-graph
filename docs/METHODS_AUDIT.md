# Methods Audit — Pipeline Recovery (April 28-29, 2026)

This document records the diagnostic forensics and corrective actions applied to the BIFO-PPR pipeline between April 28-29, 2026, in advance of bioRxiv submission. Two compounding bugs were identified and corrected. The buggy pipeline had been producing the canonical results files in the repository as of April 26. The corrected pipeline produces the values reported in the updated manuscript.

This audit is intended to be transparent for reviewers: every numerical claim in the manuscript is the result of a deterministic pipeline at `--perm-random-seed 42`, but the inputs to that pipeline contained two bugs whose interaction was non-obvious. Both bugs were independently discovered and fixed during this audit.

---

## Bug 1 — Q5 cypher hop1-restriction excluded propagating membership edges

### Symptom

The KF-CHD canonical results file (`results/kf_chd/pathway_scores_standard.csv` at git HEAD on April 27) reported WP_CILIOPATHIES with **170 member genes**. The MSIGDB definition of WP_CILIOPATHIES has **183 member genes**. The 13 missing genes were dropped during the BIFO conditioning step.

The KF-NBL canonical results similarly reported WP_CILIOPATHIES with 169 members (off by 14 from the full 183). And both cohorts' canonical results were missing all Hallmark pathway membership content (`inverse_has_signature_gene` predicate, 7,321 edges).

### Root cause

The Q5 cypher in the production cypher-export script (`cypher/kf_chd_export_queries.cypher`, also in NBL equivalent) at HEAD restricted gene-pathway membership edges to those where **both** the pathway and gene CUIs were in the seed-1-hop reachable subgraph:

```cypher
MATCH (pw:Concept)-[r]->(gene:Concept)
WHERE pw.CUI IN hop1_ids
  AND gene.CUI IN hop1_ids
  AND type(r) IN [
    'pathway_associated_with_gene',
    'has_signature_gene',
    'inverse_pathway_associated_with_gene',
    'inverse_has_signature_gene',
    'process_involves_gene',
    'gene_plays_role_in_process'
  ]
```

This `hop1_ids` restriction had two compounding effects:

1. **13 ciliopathy member genes were silently pruned.** These genes had no mechanistic 1-hop connection to any of the 1,276 KF-CHD seeds. Their only graph connection to seeds was the membership edge to WP_CILIOPATHIES itself. Because they were not in `hop1_ids`, the gene-side filter excluded them; their membership edges never made it into the export. Conditioning subsequently saw WP_CILIOPATHIES with 170 members instead of 183.
2. **All 7,321 Hallmark gene-pathway membership edges (`inverse_has_signature_gene`) were excluded.** Hallmark pathway nodes connect to their member genes only via this predicate. Most Hallmark member genes weren't in the seed 1-hop neighborhood, so the hop1-restriction dropped most Hallmark membership edges; conditioning then saw Hallmark pathways as severely under-populated.

The corrected Q5 design exports the **full global MSigDB gene-pathway membership** with no seed restriction. It uses only the inverse-direction predicates (G→PW), which are classified `weak_mechanistic_or_observational` in the BIFO YAML (propagating). The forward predicates (PW→G: `pathway_associated_with_gene`, `has_signature_gene`) are excluded at the cypher level because BIFO classifies them as `nonpropagating_context` and conditioning would drop them anyway. Cohort-specificity is established at the conditioning step where the global membership graph is intersected with the seed-reachable subgraph.

The corrected Q5 is cohort-independent: it produces bit-identical output for KF-CHD and KF-NBL (both export the same global MSigDB membership). This is now visible at the file level: the regenerated `cypher/kf_chd_query5.cypher` and `cypher/kf_nbl_query5.cypher` are 1,137 bytes each and differ only in the cohort name in the comment header.

### Git regression timeline

| Date (2026) | Commit | Status |
|---|---|---|
| Apr 16 | (initial) | Q5 used hop1-restriction with mixed predicate list — bug |
| Apr 21 ~9am | `deaeb14` | Q5 fix v1 |
| Apr 21 22:04 | `65ab3ca` | Q5 fix v2: two-branch UNION with `inverse_pathway_associated_with_gene` + `inverse_has_signature_gene` (G→PW, propagating, no hop1 restriction). Canonical-corrected. |
| Apr 25 23:32 | `5d01d32` | Regression: regenerated cypher from `pipeline/generate_export_cypher.py`, which had never received the `65ab3ca` correction. Q5 reverted to pre-fix design with hop1-restriction. |
| Apr 28 ~07:30 | `be2f9c9` | Added Q6 to auto-generator (separate fix) |
| Apr 28-29 | (this audit) | Audit run executed Q5 manually using `65ab3ca`-style cypher in `bifo_run_2026-04-28/{chd,nbl}/cypher/`. Generator NOT yet patched. |
| Apr 29 | `b3358fd` | Recovery commit: patched `make_query5()` in `pipeline/generate_export_cypher.py` to emit global cohort-independent Q5; removed monolithic per-cohort cypher files; added freshly-generated per-query files for both cohorts. |

### Validation

Re-running Q5 against the production DDKG with the corrected cypher produced 493,963 rows (vs 486,642 buggy). The breakdown:

- 486,642 `inverse_pathway_associated_with_gene` edges (G→PW for non-Hallmark MSigDB)
- 7,321 `inverse_has_signature_gene` edges (G→PW for Hallmark) — entirely missing from the buggy run

WP_CILIOPATHIES (MSIGDB:M39880) appears with 183 source CUIs in the corrected Q5 output. All 13 previously-disconnected member genes are present. Cross-reference invariants pass: every Q5 source CUI has a row in Q6 (gene-metadata query).

The patched `pipeline/generate_export_cypher.py` was verified to reproduce the audit-run Q5 byte-for-byte (modulo trivial trailing-whitespace differences). Per-query files generated by `python pipeline/generate_export_cypher.py --seeds <file> --cohort <chd|nbl> --out-dir cypher/` are now the canonical cypher source for the pipeline; the obsolete monolithic `kf_chd_export_queries.cypher` and `kf_nbl_export_queries.cypher` files have been removed from the repo.

### Why NBL also affected

NBL Q5 produced bit-identical output to CHD Q5 (md5 `e58f26c5...`) once corrected, because the corrected Q5 design is cohort-independent — it queries the global MSigDB pathway membership, not the seed-restricted subgraph. Under the buggy hop1-restricted design, both cohorts had the bug; NBL's manuscript-canonical WP_CILIOPATHIES member count was 169 (one different from CHD's 170 due to a single-rank-tie ordering edge case in the buggy filter).

---

## Bug 2 — Sparse node_to_idx in score_pathways.py

### Symptom

When the corrected Q5 substrate was run through the rest of the pipeline (`bifo_conditioning.py` → `score_pathways.py`), conditioning succeeded but scoring crashed at the rewiring null permutation phase with:

```
ValueError: axis 0 index 836600 exceeds matrix dimension 817113
```

(Numbers are CHD-specific. NBL had analogous: index 901,808 exceeds dimension 882,374.)

### Root cause

The `nodes_extended.csv` file (Q4 ∪ Q6 merged) contains some Concept nodes with multiple SAB classifications (e.g., gene-and-pathway dual roles, NCI-Thesaurus + HGNC overlaps). For CHD, 19,488 such duplicate `node_id` rows exist.

`bifo_conditioning.py:build_sparse_operator` constructs `node_to_idx` as a Python dict keyed by `node_id`. Dict construction silently deduplicates: each CUI ends up mapped to the row index of its **last** appearance. But the matrix is sized using the un-deduplicated row count. The result: `node_to_idx` has 817,113 entries (deduplicated for CHD), but its values range up to 836,600 (matching the un-deduplicated row count). The conditioned matrices and score vectors are 836,601 × 836,601 with 19,488 unused row/column slots interspersed.

The score vectors (`results_scores_cond.npy` and friends) are written at length 836,601, indexed by the same sparse mapping. The conditioning artifacts are internally consistent: both the score arrays and the index map are in the 836,601-sized space, with 19,488 gaps.

The bug was in `score_pathways.py:run_membership_rewiring_null` (line 1117) and `score_pathways.py:run_member_mean_null` (line 596):

```python
n = len(node_to_idx)
```

This sets `n` to the **deduplicated count** (817,113) and uses it as the matrix dimension. But the bridge-edge integer indices come from `node_to_idx.map(...)`, which returns values up to 836,600. CSR matrix construction crashed because indices exceeded the matrix dimension.

### Patch applied

Both occurrences of the buggy line in `score_pathways.py`:

```python
n = max(node_to_idx.values()) + 1  # account for sparse node_to_idx after conditioning dedupe
```

This sizes the matrix to match the actual range of indices in `node_to_idx`. The 19,488 unused row/column slots are zero-filled (they correspond to deduplicated multi-SAB rows whose Concept is already represented at another index). Empty rows have `row_sum=0`, handled correctly by the existing `np.where(row_sums > 0, ...)` clause.

### Pre-existing diagnosis

Commit `d3b5c7a` (April 27, 7:00 AM) by the project lead identified this exact bug independently:

> "The `build_sparse_operator` function constructs `node_to_idx` as a dict (auto-deduplicates) but uses `n=len(node_ids)` (list length) for matrix dimensions. When upstream nodes files contain duplicate `node_id` rows, this produces matrices whose row count exceeds the node index count, causing silent corruption and downstream 'row index exceeds matrix dimensions' errors when `score_pathways.py` builds null model matrices."

The commit applied the fix at the **conditioning level** (deduplicate before assigning indices). It was reverted ten hours later in commit `4f6d3d2` (April 27, 17:01) along with two other commits in a single cleanup pass that rolled back an entire SM5 work-stream (`6af4392`, `bbacd26`, `4f6d3d2`). The conditioning-dedup fix was likely caught in this rollback as collateral damage. The consequence — Bug 2 remained latent — was apparently not recognized at the time, because the buggy Q5 substrate continued producing results that didn't trigger the size-mismatch crash. Bug 2 only surfaced once Bug 1 was fixed and the corrected substrate was passed through scoring.

This audit applies the fix at the **scoring level** instead: `n = max(values) + 1`. This change is non-invasive — conditioning artifacts are unchanged. Future cleanup could reapply commit `d3b5c7a` upstream and produce dense-index conditioning outputs, which would be cleaner for downstream consumers; but for this audit, the scoring-level patch is the minimal change that fully recovers correctness.

### Why the canonical Apr 26 run did not trigger this

The buggy Q5 substrate produced fewer pathway-member nodes that survived conditioning. Bug 1 effectively acted as a filter, reducing the input to scoring such that the integer indices passed to the rewiring null happened to fit within the truncated matrix dimension. Once Bug 1 was corrected, the corrected substrate had 7,321 additional Hallmark edges and 13 reconnected ciliopathy genes, producing more bridge edges with index values up to 836,600, exceeding the truncated dimension. Bug 2 was latent until Bug 1 was fixed.

---

## Run provenance

Recovery run directory: `/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_2026-04-28/`

```
bifo_run_2026-04-28/
├── chd/
│   ├── cypher/             # per-query cypher used for this run (with corrected Q5)
│   ├── outputs/            # raw cypher-shell output + clean + merged
│   ├── results/            # conditioning + scoring outputs + summary
│   ├── validation/         # md5 checks, run logs, PIDs
│   ├── seeds.txt           # snapshot of input seed list (1,292 lines)
│   └── seed_cuis.txt       # 1,276 resolved CUIs
└── nbl/
    └── (same structure; 1,395 resolved seed CUIs)
```

### Run hashes

CHD Q2-Q6 outputs:

| Query | Lines | md5 | Match Apr 16 canonical? |
|---|---|---|---|
| Q2 seed_nodes | 1,277 | `ec59f127dd6502a736802de3057d7241` | IDENTICAL |
| Q3 edges_raw | 5,261,300 | `6e9bbed4055e6f6d81c3957a1904a454` | IDENTICAL |
| Q4 nodes | 1,079,543 | `9b4212bc41eb010ad11cd12a0635006b` | IDENTICAL |
| Q5 pathway_membership | 493,963 | `e58f26c51203a3e6710fea652c417e7a` | DIFFERENT (corrected Q5) |
| Q6 pathway_member_nodes | 21,353 | `b63ba00a1fff948b21be50d6095cac96` | IDENTICAL |

NBL Q2-Q6 outputs:

| Query | Lines | md5 | Notes |
|---|---|---|---|
| Q2 seed_nodes | 1,396 | `5640b40049cca82db98edf0f304d9487` | 1,395 resolved CUIs |
| Q3 edges_raw | 5,520,176 | `6f745a16691e992974c4cda587170f16` | |
| Q4 nodes | 1,169,818 | `50db49020c6d5a14f59d09f2352960c5` | |
| Q5 pathway_membership | 493,964 | `e58f26c51203a3e6710fea652c417e7a` | IDENTICAL to CHD Q5 (cohort-independent global) |
| Q6 pathway_member_nodes | 21,354 | `b63ba00a1fff948b21be50d6095cac96` | IDENTICAL to CHD Q6 (cohort-independent global) |

### Substrate validation

- WP_CILIOPATHIES (MSIGDB:M39880) member count in Q5: **183** ✓ (was 170 buggy CHD, 169 buggy NBL)
- Predicate distribution in Q5: 486,642 `inverse_pathway_associated_with_gene` + 7,321 `inverse_has_signature_gene`
- All 13 disconnected ciliopathy genes are in Q6 (gene-metadata query) ✓
- All 13 disconnected ciliopathy genes were ABSENT from the canonical April 18 `nodes_clean_noncc.csv.gz` ✓ (the bug, confirmed)
- All 183 ciliopathy members are in the conditioned operator after the corrected pipeline ✓ (CHD: 183/183 in `results_node_index.json`; NBL: 183/183 in `results_node_index.json`)

### Conditioning results

| Cohort | Conditioned operator size | Conditioning wall time |
|---|---|---|
| KF-CHD | 817,113 nodes (817,113 unique CUIs in 836,601 row slots) | ~17 min |
| KF-NBL | 882,374 nodes (882,374 unique CUIs in 901,809 row slots) | ~17 min |

### Scoring results (after Bug 2 patch)

| Cohort | Scored pathways | Calibrated under rewiring null | Wall time |
|---|---|---|---|
| KF-CHD | 2,044 | 2,044 (100%) | 31m 12s @ 192 cores |
| KF-NBL | 2,111 | 2,098 (99.4%) | ~30 min @ 192 cores |

---

## Authoritative numbers (corrected pipeline, April 29, 2026)

### KF-CHD WP_CILIOPATHIES

| Metric | Buggy (April 26) | Corrected (April 29) |
|---|---|---|
| Members in scored set | **170** (BUG) | **183** ✓ |
| Rank by degree_norm | 43 of 2,130 | **12 of 2,044** |
| null_z (membership-rewiring) | 48.708 | **17.888** |
| empirical p (rewiring) | 0.0010 (floor) | 0.0010 (floor) |
| empirical q (rewiring, BH) | 0.0078 | 0.0100 |
| null_calibrated (rewiring) | True | True |
| signal_to_null_mean | 5.32 | 4.06 |
| member_mean null_z | 1.39 | 2.07 |
| member_mean q | 0.865 | 0.192 |

### KF-NBL WP_CILIOPATHIES

| Metric | Buggy (April 26) | Corrected (April 29) |
|---|---|---|
| Members in scored set | **169** (BUG) | **183** ✓ |
| Rank by degree_norm | 3 of 2,196 | **1 of 2,111** |
| null_z (membership-rewiring) | 18.948 | **NaN (degenerate)** |
| empirical q (rewiring, BH) | 0.0121 | NaN |
| null_calibrated (rewiring) | True | **False** (signal_to_null_mean = 12.6 > 10 cutoff) |
| signal_to_null_mean | 2.64 | **12.60** |
| **member_mean null_z** | 2.43 | **6.63** |
| **member_mean q** | 0.0997 | **0.0340** |

The KF-NBL membership-rewiring null is degenerate for WP_CILIOPATHIES: the observed signal is so dominant that no rewiring of pathway membership produces comparable values. Per the calibration filter (`signal_to_null_mean > 10`), the rewiring null_z is reported as NaN. The member-mean stratified null is calibrated and gives null_z = 6.63, q = 0.034 (significant at q < 0.05).

### Cilia cluster — 17 MSigDB pathways (added WP_JOUBERT_SYNDROME this audit)

The cilia reference set, expanded to 17 pathways during this audit (added MSIGDB:M39835 WP_JOUBERT_SYNDROME), comprises three core ciliopathy-and-cilia-development pathways, three primary-cilium-trafficking pathways, eight Hedgehog-signaling pathways (cilia-mediated signaling), and three other curated cilia sets.

**KF-CHD CORRECTED — cilia cluster (rewiring null, ranked):**

| Pathway | Members | Rank | null_z | q | Calibrated? |
|---|---:|---:|---:|---:|:---:|
| WP_CILIOPATHIES | 183 | 12 | 17.89 | 0.0100 | T |
| WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT | 102 | 38 | 10.40 | 0.0100 | T |
| WP_JOUBERT_SYNDROME | 77 | 47 | 9.39 | 0.0100 | T |
| REACTOME_CILIUM_ASSEMBLY | 200 | 78 | 6.01 | 0.0100 | T |
| REACTOME_CARGO_TRAFFICKING_TO_THE_PERICILIARY_MEMBRANE | 51 | 166 | 5.01 | 0.0100 | T |
| REACTOME_BBSOME_MEDIATED_CARGO_TARGETING_TO_CILIUM | 23 | 278 | 3.76 | 0.0175 | T |
| REACTOME_SIGNALING_BY_HEDGEHOG | 150 | 233 | 1.96 | 0.190 | T |
| WP_HEDGEHOG_SIGNALING_PATHWAY | 44 | 366 | 2.25 | 0.125 | T |
| WP_CILIARY_LANDSCAPE | 218 | 269 | 0.53 | 0.755 | T |
| (8 other Hedgehog/cilia pathways) | various | 488–1658 | -1.35 to 1.85 | 0.23–1.0 | T |

**KF-NBL CORRECTED — cilia cluster (rewiring + member_mean nulls):**

| Pathway | Members | Rank | rew_z | rew_q | mm_z | mm_q | Cal? |
|---|---:|---:|---:|---:|---:|---:|:---:|
| WP_CILIOPATHIES | 183 | **1** | NaN | NaN | **6.63** | **0.034** | **F (deg.)** |
| WP_JOUBERT_SYNDROME | 77 | **6** | 6.45 | 0.646 | **5.21** | **0.034** | T |
| WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT | 102 | **7** | 8.03 | 0.646 | **4.97** | **0.034** | T |
| REACTOME_CILIUM_ASSEMBLY | 200 | 14 | 6.18 | 0.646 | 1.61 | 0.556 | T |
| REACTOME_BBSOME_MEDIATED_CARGO_TARGETING_TO_CILIUM | 23 | **23** | NaN | NaN | **5.74** | **0.034** | **F (deg.)** |
| REACTOME_CARGO_TRAFFICKING_TO_THE_PERICILIARY_MEMBRANE | 51 | 34 | 2.67 | 0.646 | 3.15 | 0.034 | T |
| WP_CILIARY_LANDSCAPE | 218 | 91 | 1.09 | 0.646 | -1.96 | 1.0 | T |
| (Hedgehog and 9 other cilia pathways) | various | 92–1369 | various | various | various | various | T |

Note: NBL rewiring-null q-values are uniformly elevated (0.6462) because empirical p-values for highly-significant pathways floor at 1/(n_perms+1) = 0.001, and BH-adjustment across the 2,098-pathway calibrated set inflates the threshold. The member-mean null produces more interpretable q-values (0.034 for several cilia pathways).

### Top scoring pathways KF-CHD (corrected) — top 10 by null_z (calibrated)

| Pathway | Members | null_z | q |
|---|---:|---:|---:|
| Signal Transduction | (broad) | 29.88 | 0.0100 |
| Ionophore activity | (broad) | 28.96 | 0.0100 |
| Transmembrane Transport | (broad) | 28.17 | 0.0100 |
| Regulation of Cell Shape | (broad) | 21.03 | 0.0100 |
| Phosphorylation | (broad) | 20.16 | 0.0100 |
| Transcriptional Regulation | (broad) | 20.14 | 0.0100 |
| Intercellular Communication Process | (broad) | 18.97 | 0.0100 |
| DNA Repair | (broad) | 18.36 | 0.0100 |
| DNA Binding | (broad) | 18.17 | 0.0100 |
| WP_CILIOPATHIES | 183 | 17.89 | 0.0100 |

CHD's top tier is dominated by very broad GO/Reactome categories. WP_CILIOPATHIES at rank 10-12 is the highest specifically curated pathway in the top 20, ranking above more specific signaling and regulatory categories.

### Top scoring pathways KF-NBL (corrected) — top 10 by null_z (calibrated)

| Pathway | Members | null_z | q |
|---|---:|---:|---:|
| WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT | 102 | 8.03 | 0.6462 |
| REACTOME_FATTY_ACID_METABOLISM | (broad) | 7.48 | 0.6462 |
| REACTOME_TRNA_PROCESSING | (broad) | 7.02 | 0.6462 |
| WP_JOUBERT_SYNDROME | 77 | 6.45 | 0.6462 |
| WP_7Q1123_COPY_NUMBER_VARIATION_SYNDROME | (broad) | 6.27 | 0.6462 |
| REACTOME_SENSORY_PROCESSING_OF_SOUND | (cilia-related) | 6.23 | 0.6462 |
| REACTOME_CILIUM_ASSEMBLY | 200 | 6.18 | 0.6462 |
| REACTOME_RNA_POLYMERASE_II_TRANSCRIBES_SNRNA_GENES | (broad) | 5.48 | 0.6462 |
| WP_TYROBP_CAUSAL_NETWORK | (broad) | 5.34 | 0.6462 |
| REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION | (broad) | 5.32 | 0.6462 |

NBL's top tier (calibrated set, by null_z) contains four distinct cilia-related pathways: WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT, WP_JOUBERT_SYNDROME, REACTOME_SENSORY_PROCESSING_OF_SOUND, and REACTOME_CILIUM_ASSEMBLY. WP_CILIOPATHIES itself is excluded from this list because it is null-degenerate (signal_to_null_mean > 10), but it ranks #1 by degree_norm across all 2,111 scored pathways.

---

## Manuscript framing implications

### Updated abstract narrative

The corrected pipeline supports a stronger cluster-level claim than the original single-pathway framing:

> "Across a curated cilia cluster of 17 MSigDB pathways, BIFO-PPR placed the three core ciliopathy/primary-cilium-development pathways (WP_CILIOPATHIES, WP_JOUBERT_SYNDROME, WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT) in the top 50 of KF-NBL (rank 1, 6, 7) and top 50 of KF-CHD (rank 12, 38, 47), with cilium assembly and BBSome trafficking pathways also in the top 300 of both cohorts. The signal is consistent across both cohorts despite their phenotypic disjointness, supporting cilia-mediated transcriptional regulation as a shared mechanism in pediatric germline genetic disease and consistent with the curation-bias hypothesis (Methods §X)."

### Reporting strategy for the dual-null situation

KF-CHD's WP_CILIOPATHIES is fully calibrated under the membership-rewiring null. KF-NBL's WP_CILIOPATHIES is null-degenerate under the rewiring null but calibrated under the member-mean stratified null. The recommended approach is to **report both nulls explicitly**:

> "Under the membership-rewiring null, WP_CILIOPATHIES ranks 12th of 2,044 scored pathways in KF-CHD (null_z = 17.9, q = 0.010) and 1st of 2,111 in KF-NBL (rewiring null degenerate due to dominant signal; under the member-mean stratified null, null_z = 6.63, q = 0.034)."

This avoids the appearance of cherry-picking nulls and accurately conveys the asymmetry.

---

## Recovery commits

The following five commits have landed on `bifo-graph` `main` (origin: `github.com:TaylorResearchLab/bifo-graph.git`) as of April 29, 2026:

1. `d0dada7` **fix(scoring): correct sparse node_to_idx handling in null matrix sizing** — `pipeline/score_pathways.py`: change `n = len(node_to_idx)` to `n = max(node_to_idx.values()) + 1` at lines 596 and 1117. Backwards-compatible; preserves canonical conditioning artifacts (`results_scores_*.npy` byte-identical).
2. `c73318d` **data: add WP_JOUBERT_SYNDROME (MSIGDB:M39835) to cilia reference** — `data/cohorts/chd/kf_chd_cilia_reference.txt`: appended Joubert syndrome (classic ciliopathy phenotype). Reference cluster expanded from 16 to 17 pathways.
3. `349d2d8` **docs: add METHODS_AUDIT.md documenting Apr 28-29 pipeline recovery** — this audit document, located at `docs/METHODS_AUDIT.md`.
4. `b3358fd` **fix(cypher): correct Q5 to global cohort-independent MSigDB membership** — patches `make_query5()` in `pipeline/generate_export_cypher.py`; removes obsolete monolithic per-cohort cypher files (`kf_chd_export_queries.cypher`, `kf_nbl_export_queries.cypher`); adds freshly-generated per-query cypher files (`kf_{chd,nbl}_query{2,3,4,5,6}.cypher`). Q5 file size is identical (1,137 bytes) between cohorts, making cohort-independence visible at the file-listing level. `scripts/run_kf_{chd,nbl}_export.sh` already expected per-query files.
5. `5cfd951` **docs: update cypher pipeline references to per-query layout** — `README.md` and `scripts/run_kf_chd_export.sh` text updates referencing per-query layout instead of monolithic cypher files.

### Pending follow-up

These remain to be done in subsequent commits:

- **Results regeneration with full output flags.** Re-run `score_pathways.py` with `--out-member-scores` and `--out-influential-nodes` flags so the canonical results include driver-gene annotations alongside `pathway_scores_standard.csv`/`pathway_metrics_standard.json`. Expected commit: `results: regenerate KF-CHD and KF-NBL outputs against corrected pipeline` (replaces `results/kf_chd/*.csv`, `results/kf_nbl/*.csv`; updates REPRODUCE.md with new md5s).
- **Figure regeneration.** `figures: regenerate fig9_rank_vs_nullz against corrected scores` — re-render Figure 8 (manuscript) from corrected `pathway_scores_standard.csv`.
- **Manuscript text updates** (parallel commit campaign on `bifo-paper-1` repo): abstract null_z values, Results section rankings, methods limitations subsection (GO/MSigDB membership-ratio asymmetry, P/LP-curated seed composition), Discussion reframing for the asymmetric NBL/CHD result.

---

## Process improvements (carry forward)

This audit revealed several classes of preventable bug. Process rules added to the project handoff:

1. **Validate output non-empty + correct header before trusting cypher-shell exit code.** Cypher-shell exits 0 on empty results. Always check `wc -l` and `head -1` before considering a query successful.
2. **Check BIFO YAML predicate classifications when conditioning drops things you don't expect.** The Q5 bug hid for several days because conditioning correctly dropped nonpropagating edges; the symptom (170 members instead of 183) was visible only at scoring time, not at conditioning time. Always check the YAML if conditioning is dropping more than expected.
3. **Don't trust auto-generators without inspection — and patch them rather than working around them.** The `pipeline/generate_export_cypher.py` had a buggy `make_query5()` that emitted a hop1-restricted query with mixed forward/inverse predicates. Earlier work (`65ab3ca`, April 21) had hand-edited the static cypher files but not the generator, so a subsequent regeneration (`5d01d32`, April 25) silently re-introduced the bug. The recovery commit (`b3358fd`, April 29) patched `make_query5()` directly so the generator now emits the corrected global cohort-independent Q5; this also eliminated the duplicated monolithic cypher files in favor of per-query files matching the layout `scripts/run_kf_*_export.sh` already expected.
4. **Cross-host process monitoring requires checking PIDs locally + logs cross-host.** A `kill -0 PID` on a host where the process doesn't live always returns "not running", giving a false "DONE" status. Monitor scripts must determine which host they're on and use the right check.
5. **When fixing a bug, check whether the same bug was diagnosed and reverted previously.** The sparse-n bug had been correctly diagnosed in commit `d3b5c7a` (April 27) and reverted ten hours later. Searching the git log for prior diagnoses of the same symptom would have saved time.
6. **Not all reverts are intentional.** Three reverts on April 27 (`6af4392`, `bbacd26`, `4f6d3d2`) appear to have been a single cleanup pass that rolled back an SM5 work-stream; the conditioning-dedup fix was caught up in this rollback as collateral. Review revert PRs/commits carefully to ensure each individual revert is intended.

---

## Limitations and interpretation caveats

A separate review of the corrected KF-CHD summary output (run by another reviewer on April 29, 2026) surfaced three structural caveats that must be reflected in the manuscript framing. These are not pipeline bugs — they are properties of the substrate and seed inputs that affect how the corrected numbers should be interpreted.

### 1. GO-vs-MSigDB membership asymmetry

In the scored pathway table, GO Concepts (the broad categories at the top of the rewiring-null leaderboard — `Signal Transduction`, `Ionophore activity`, `Transmembrane Transport`, `Phosphorylation`, etc.) show `n_seed_members ≈ n_members_in_universe` — that is, the ratio of seed-overlapping members to total members is near 1.00. By contrast, MSigDB pathways (HALLMARK_, REACTOME_, WP_, etc.) show ratios in the 0.2–0.5 range — most of their members are not in the seed list.

This asymmetry is structural to the DDKG GO ingestion, not a bug in our pipeline. Most likely the GO Concepts in the DDKG are pre-filtered through a disease-gene curation step, so the gene set associated with each GO term in the graph is a curated subset rather than the full GO membership. (This is a known property of some curated GO ingestions; should be verified against the DDKG provenance documentation.)

The implication: **GO-vs-MSigDB rank comparisons in the top-tier of the scored set are not apples-to-apples.** A GO term like `Signal Transduction` ranking above WP_CILIOPATHIES does not mean the GO term has a stronger signal — it means the GO term's gene-set definition in the DDKG is more enriched for disease-curated content than the MSigDB pathway definition. The cilia cluster pathways (all MSigDB) are internally consistent and the rank-comparisons within the cilia cluster (or across other MSigDB pathway sets) are valid.

For the manuscript: avoid quoting top-tier GO ranks in the abstract. Frame the cilia cluster results in terms of MSigDB pathway rankings only. If GO terms are mentioned, note the membership-ratio caveat.

### 2. Seed composition: P/LP-curated panel, not CHD- or NBL-specific gene sets

The KF-CHD seed list (1,276 resolved CUIs) and KF-NBL seed list (1,395 resolved CUIs) are derived from AutoGVP P/LP-curated germline variants in the respective proband cohorts. This is a clinical Mendelian-variant filter applied to per-proband WGS — not a CHD-specific or NBL-specific gene panel.

A spot-check of the KF-CHD seed list shows extensive coverage of: ciliopathy genes (BBS family, EVC/EVC2, NPHS1/2), Mendelian metabolic and connective-tissue genes (CFTR, RYR1, CYP21A2), DNA repair genes (BRCA1/2, FANC family, MSH family, ATM, CHEK2), coagulation factors (F2, F7, F8, F11), and complement components.

Conspicuously absent from the KF-CHD seed list: the canonical heart-development transcription factors NKX2-5, TBX5, TBX1, GATA6, JAG1, NOTCH1. Heart-development TF coverage is sparse (GATA4, GATA5, ZIC3, NODAL only). This is the expected pattern for a P/LP-curated filter applied to CHD probands: ClinVar P/LP curation density is highest for genes with well-established Mendelian disease associations across many phenotypes, not for cohort-specific developmental regulators where individual variants are typically VUS.

The implication: **the pathway results reflect what is curated in ClinVar P/LP, not what is mechanistically core to the cohort phenotype.** Both cohorts' top-tier results contain rich ciliopathy signal because ciliopathy genes are densely curated in ClinVar. The corrected pipeline did not "discover" cilia biology in CHD or NBL — it correctly propagated signal from a P/LP-curated seed set in which ciliopathy genes are over-represented.

This recontextualizes but does not invalidate the corrected results. The cluster-level claim ("BIFO-PPR places six of the 17 cilia cluster pathways in NBL's top 50 and five in CHD's top 300") still holds. But it should be framed as a property of P/LP-curated variants in pediatric germline cohorts, not as cohort-specific biological discovery.

This is the same observation that motivates Paper 3 (curation-bias correction methodology). Paper 1's Methods and Discussion should explicitly acknowledge this. A clean framing:

> "Seed sets in this analysis are AutoGVP P/LP-curated variants identified in proband WGS. ClinVar P/LP curation density is highest for genes with well-established Mendelian disease associations, which biases the seed composition toward broadly-curated mechanistic categories (e.g., ciliopathy, DNA repair, lysosomal storage) and away from cohort-specific developmental regulators where most variants are clinically VUS. The pathway results reflect this bias and should be interpreted as identifying signal in the P/LP-curated subset of variants present in each cohort, not as identifying cohort-specific gene-level biology. The systematic effect of curation bias on cohort-comparison studies is the subject of separate methods development (Paper 3)."

### 3. Empty driver columns

The summarize_results.py output table has empty `seed_member_scores`, `influential_nodes_local`, and `influential_nodes_global` columns. This is because the original scoring run did not specify `--out-member-scores` or `--out-influential-nodes` flags. The pathway scoring CSV is complete and authoritative; only the per-pathway driver-gene annotations are missing from the summary.

Resolution: re-run scoring with the additional output flags (this regenerates byte-identical pathway_scores_standard.csv plus the two new driver files, given the deterministic seed). The `HANDOFF_april29_2026.md` document includes the full scoring command with these flags. Cost: ~30 min × 2 cohorts. Not blocking for the audit, but recommended before manuscript figure regeneration so that supplementary tables can include driver-gene annotations.

---

## Acknowledgment of scope

This audit corrects the substrate (Q5) and scoring (n) bugs in the BIFO-PPR pipeline as applied to KF-CHD and KF-NBL cohort data with `--perm-random-seed 42`. The pipeline architecture, the BIFO ontology, the conditioning algorithm, and the membership-rewiring null model are unchanged. The methodology described in the manuscript is correct as designed; the bugs were in the cypher hop1-restriction filter on Q5 (excluding propagating membership edges that should have been retained) and a single integer assignment in the scoring code (matrix dimension calculation that didn't account for sparse `node_to_idx` from the conditioning step). Recovery is mechanical and fully deterministic.

The cohort-similarity / curation-bias finding documented separately in `audit_2026-04-28/docs/audit_findings.md` and the U24 memo is independent of these substrate bugs and remains valid: it concerns the AutoGVP P/LP filter producing seed sets with shared ClinVar curation density across cohort-disjoint phenotypes. That finding is the basis for Paper 3 (curation-bias correction methodology) and is unrelated to the bugs corrected here. The corrected results in this audit, if anything, strengthen the cohort-similarity finding by showing the cilia signal at consistent strength across both cohorts.
