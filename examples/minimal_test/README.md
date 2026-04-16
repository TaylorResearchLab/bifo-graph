# Minimal BIFO Pipeline Test Run

End-to-end validation with a small CHD seed set.
Goal: confirm the pipeline runs and produces non-trivial output —
biological correctness comes later with the full cohort gene list.

## Files in this directory

```
neo4j_export.cypher        — Export query: 2-hop CHD neighborhood
seed_nodes.txt             — 10 known CHD genes as seeds (HGNC CUIs)
heldout_nodes.txt          — 5 held-out CHD genes for evaluation
run_test.sh                — One-command pipeline run
```

## Step 1 — Export from Neo4j

Run `neo4j_export.cypher` against your DDKG instance (Neo4j Browser or
cypher-shell). It exports two CSV files:
  - nodes.csv
  - edges_raw.csv

See the query file for exact output format requirements.

## Step 2 — Run the pipeline

```bash
cd minimal_test_run/
bash run_test.sh
```

This runs bifo_conditioning.py and then score_pathways.py.
Outputs written to ./test_output/:
  - test_output/results.json          — four-arm PPR metrics (gene recovery)
  - test_output/results_scores_cond.npy
  - test_output/results_scores_raw.npy
  - test_output/results_node_index.json
  - test_output/results_kept_edges.csv
  - test_output/pathway_scores.csv    — Analysis B ranked pathways
  - test_output/pathway_metrics.json

## What success looks like

- results.json has non-NaN auroc/auprc values
- pathway_scores.csv has rows (pathway nodes were found)
- raw_vs_cond_delta is non-zero for top pathways (conditioning has effect)
- rank_improvement_cond_vs_raw is positive (conditioning helps)

## Seed and heldout gene list

Seeds and heldout are expressed as HGNC Concept node IDs (CUIs).
The query below retrieves the correct CUI for any HGNC gene symbol:

  MATCH (c:Concept)-[:CODE]->(code:Code)
  WHERE code.SAB = 'HGNC' AND code.CodeID = 'HGNC:XXXX'
  RETURN c.CUI AS cui

Or by symbol via PREF_TERM:
  MATCH (c:Concept)-[:PREF_TERM]->(t)
  MATCH (c)-[:CODE]->(code:Code)
  WHERE code.SAB = 'HGNC' AND t.name = 'GATA4'
  RETURN c.CUI AS cui LIMIT 1
