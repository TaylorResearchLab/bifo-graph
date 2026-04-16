// =============================================================================
// BIFO C4 Control Export — Notch Signaling Control
// Seeds: 30 CUIs  |  Heldout: 14 CUIs  |  Total: 44 CUIs
//
// These are the C4 randomization control analyses for the CHD curated benchmark.
// The graph is exported using all CUIs (seeds + heldout) as the neighborhood center,
// then bifo_conditioning.py is run with --seed-nodes and --heldout-nodes separately.
//
// Reference file: data/benchmark/c4_notch_pathway_reference.txt
//
// Output columns match bifo_conditioning.py:
//   nodes.csv     — node_id, label, name, sab
//   edges_raw.csv — source, target, predicate
// =============================================================================

// ---------------------------------------------------------------------------
// QUERY 1 — Export 1-hop nodes → save as: c4_notch_nodes.csv
// Columns: node_id, label, name, sab
// Uses CASE ordering to prefer HGNC SAB over other SABs
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C0812258', 'C1412637', 'C1336607', 'C1415524', 'C1336558', 'C1421267',
  'C1415526', 'C1417767', 'C1415725', 'C1414070', 'C1412638', 'C1367721',
  'C1418285', 'C1334889', 'C1415525', 'C1334498', 'C1423876', 'C1538913',
  'C0812223', 'C1419030', 'C1412137', 'C1424603', 'C1366763', 'C1416828',
  'C1417768', 'C1424142', 'C1335823', 'C1415380', 'C1417230', 'C1414072',
  'C1412719', 'C1416171', 'C1419044', 'C0812228', 'C1416525', 'C1333569',
  'C1414224', 'C1419041', 'C1418921', 'C1537839', 'C1419037', 'C1823863',
  'C1823300', 'C1414069'
]
MATCH (seed)-[r1]-(hop1:Concept)
WHERE type(r1) <> 'CODE'
  AND type(r1) <> 'PREF_TERM'
  AND type(r1) <> 'STY'

WITH collect(DISTINCT seed) + collect(DISTINCT hop1) AS all_concepts

UNWIND all_concepts AS c
WITH DISTINCT c

OPTIONAL MATCH (c)-[:CODE]->(code:Code)
WITH c, code
ORDER BY
  CASE code.SAB
    WHEN 'HGNC'     THEN 0
    WHEN 'NCBIGENE' THEN 1
    WHEN 'MSIGDB'   THEN 2
    WHEN 'REACTOME' THEN 3
    WHEN 'GO'       THEN 4
    WHEN 'MONDO'    THEN 5
    WHEN 'OMIM'     THEN 6
    ELSE 99
  END ASC
WITH c, collect(code)[0] AS best_code
OPTIONAL MATCH (c)-[:PREF_TERM]->(pt:Term)

RETURN
  c.CUI                              AS node_id,
  'Concept'                          AS label,
  coalesce(pt.name, c.CUI)           AS name,
  coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;


// ---------------------------------------------------------------------------
// QUERY 2 — Export 1-hop edges → save as: c4_notch_edges_raw.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C0812258', 'C1412637', 'C1336607', 'C1415524', 'C1336558', 'C1421267',
  'C1415526', 'C1417767', 'C1415725', 'C1414070', 'C1412638', 'C1367721',
  'C1418285', 'C1334889', 'C1415525', 'C1334498', 'C1423876', 'C1538913',
  'C0812223', 'C1419030', 'C1412137', 'C1424603', 'C1366763', 'C1416828',
  'C1417768', 'C1424142', 'C1335823', 'C1415380', 'C1417230', 'C1414072',
  'C1412719', 'C1416171', 'C1419044', 'C0812228', 'C1416525', 'C1333569',
  'C1414224', 'C1419041', 'C1418921', 'C1537839', 'C1419037', 'C1823863',
  'C1823300', 'C1414069'
]

MATCH (seed)-[r]-(neighbor:Concept)
WHERE type(r) <> 'CODE'
  AND type(r) <> 'PREF_TERM'
  AND type(r) <> 'STY'

RETURN DISTINCT
  seed.CUI     AS source,
  neighbor.CUI AS target,
  type(r)      AS predicate
ORDER BY source, predicate, target;


// ---------------------------------------------------------------------------
// QUERY 3 — Pathway membership edges → save as: c4_notch_pathway_membership_edges.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C0812258', 'C1412637', 'C1336607', 'C1415524', 'C1336558', 'C1421267',
  'C1415526', 'C1417767', 'C1415725', 'C1414070', 'C1412638', 'C1367721',
  'C1418285', 'C1334889', 'C1415525', 'C1334498', 'C1423876', 'C1538913',
  'C0812223', 'C1419030', 'C1412137', 'C1424603', 'C1366763', 'C1416828',
  'C1417768', 'C1424142', 'C1335823', 'C1415380', 'C1417230', 'C1414072',
  'C1412719', 'C1416171', 'C1419044', 'C0812228', 'C1416525', 'C1333569',
  'C1414224', 'C1419041', 'C1418921', 'C1537839', 'C1419037', 'C1823863',
  'C1823300', 'C1414069'
]

MATCH (seed)-[r1]-(neighbor:Concept)
WHERE type(r1) <> 'CODE'
  AND type(r1) <> 'PREF_TERM'
  AND type(r1) <> 'STY'

WITH collect(DISTINCT seed) + collect(DISTINCT neighbor) AS all_genes

UNWIND all_genes AS gene
MATCH (gene)-[r]-(pathway:Concept)
WHERE type(r) IN [
    'pathway_associated_with_gene',
    'inverse_pathway_associated_with_gene',
    'has_signature_gene',
    'inverse_has_signature_gene',
    'process_involves_gene',
    'gene_plays_role_in_process'
]
RETURN DISTINCT
    gene.CUI    AS source,
    pathway.CUI AS target,
    type(r)     AS predicate
ORDER BY source, predicate, target;


// =============================================================================
// RUN COMMANDS (after export + clean_cypher_output.py on each file):
//
//   python bifo_conditioning.py \
//     --nodes         c4_notch_nodes.csv \
//     --edges         c4_notch_edges_all.csv \
//     --mapping       config/bifo_ddkg_mapping.yaml \
//     --seed-nodes    data/benchmark/c4_notch_seed_nodes.txt \
//     --heldout-nodes data/benchmark/c4_notch_heldout_nodes.txt \
//     --out-json      results/c4_notch/results.json
//
//   python score_pathways.py \
//     --nodes             c4_notch_nodes.csv \
//     --edges-raw         c4_notch_edges_all.csv \
//     --edges-conditioned results/c4_notch/results_kept_edges.csv \
//     --scores-cond       results/c4_notch/results_scores_cond.npy \
//     --scores-raw        results/c4_notch/results_scores_raw.npy \
//     --node-index        results/c4_notch/results_node_index.json \
//     --seed-nodes        data/benchmark/c4_notch_seed_nodes.txt \
//     --chd-pathways      data/benchmark/c4_notch_pathway_reference.txt \
//     --out-csv           results/c4_notch/pathway_scores.csv \
//     --out-json          results/c4_notch/pathway_metrics.json
// =============================================================================
