// =============================================================================
// BIFO C4 Control Export — MAPK/RAS Signaling Control
// Seeds: 63 CUIs  |  Heldout: 28 CUIs  |  Total: 91 CUIs
//
// These are the C4 randomization control analyses for the CHD curated benchmark.
// The graph is exported using all CUIs (seeds + heldout) as the neighborhood center,
// then bifo_conditioning.py is run with --seed-nodes and --heldout-nodes separately.
//
// Reference file: data/benchmark/c4_mapk_pathway_reference.txt
//
// Output columns match bifo_conditioning.py:
//   nodes.csv     — node_id, label, name, sab
//   edges_raw.csv — source, target, predicate
// =============================================================================

// ---------------------------------------------------------------------------
// QUERY 1 — Export 1-hop nodes → save as: c4_mapk_nodes.csv
// Columns: node_id, label, name, sab
// Uses CASE ordering to prefer HGNC SAB over other SABs
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1538917', 'C0079471', 'C1418844', 'C1415791', 'C0812296', 'C1422215',
  'C1825576', 'C1335280', 'C1332109', 'C1333579', 'C1539630', 'C0919506',
  'C1414188', 'C1825057', 'C1335821', 'C1334290', 'C1419279', 'C1325331',
  'C1332721', 'C0879590', 'C1334291', 'C1333539', 'C1367721', 'C1334479',
  'C0919466', 'C1335284', 'C1418901', 'C0812215', 'C1970159', 'C1366644',
  'C1333572', 'C1333669', 'C1419044', 'C1420314', 'C1419030', 'C1335212',
  'C1334140', 'C0812265', 'C1332790', 'C1333535', 'C1425013', 'C1334522',
  'C0919507', 'C1414607', 'C1412166', 'C1418926', 'C0080298', 'C1420385',
  'C1414603', 'C1418576', 'C1334117', 'C1336620', 'C1414593', 'C1414592',
  'C1336929', 'C1538740', 'C1335202', 'C1334906', 'C1422166', 'C1332791',
  'C1335201', 'C1334480', 'C1417230', 'C1334863', 'C1425540', 'C1367472',
  'C1366894', 'C1425048', 'C0242957', 'C1335823', 'C1416498', 'C1421267',
  'C1419280', 'C1332453', 'C1333707', 'C1414655', 'C1333928', 'C1418903',
  'C1333578', 'C1366765', 'C1334292', 'C1857782', 'C1333362', 'C1419282',
  'C1417668', 'C1419037', 'C1419041', 'C1333543', 'C1367476', 'C0812290',
  'C1334116'
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
// QUERY 2 — Export 1-hop edges → save as: c4_mapk_edges_raw.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1538917', 'C0079471', 'C1418844', 'C1415791', 'C0812296', 'C1422215',
  'C1825576', 'C1335280', 'C1332109', 'C1333579', 'C1539630', 'C0919506',
  'C1414188', 'C1825057', 'C1335821', 'C1334290', 'C1419279', 'C1325331',
  'C1332721', 'C0879590', 'C1334291', 'C1333539', 'C1367721', 'C1334479',
  'C0919466', 'C1335284', 'C1418901', 'C0812215', 'C1970159', 'C1366644',
  'C1333572', 'C1333669', 'C1419044', 'C1420314', 'C1419030', 'C1335212',
  'C1334140', 'C0812265', 'C1332790', 'C1333535', 'C1425013', 'C1334522',
  'C0919507', 'C1414607', 'C1412166', 'C1418926', 'C0080298', 'C1420385',
  'C1414603', 'C1418576', 'C1334117', 'C1336620', 'C1414593', 'C1414592',
  'C1336929', 'C1538740', 'C1335202', 'C1334906', 'C1422166', 'C1332791',
  'C1335201', 'C1334480', 'C1417230', 'C1334863', 'C1425540', 'C1367472',
  'C1366894', 'C1425048', 'C0242957', 'C1335823', 'C1416498', 'C1421267',
  'C1419280', 'C1332453', 'C1333707', 'C1414655', 'C1333928', 'C1418903',
  'C1333578', 'C1366765', 'C1334292', 'C1857782', 'C1333362', 'C1419282',
  'C1417668', 'C1419037', 'C1419041', 'C1333543', 'C1367476', 'C0812290',
  'C1334116'
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
// QUERY 3 — Pathway membership edges → save as: c4_mapk_pathway_membership_edges.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1538917', 'C0079471', 'C1418844', 'C1415791', 'C0812296', 'C1422215',
  'C1825576', 'C1335280', 'C1332109', 'C1333579', 'C1539630', 'C0919506',
  'C1414188', 'C1825057', 'C1335821', 'C1334290', 'C1419279', 'C1325331',
  'C1332721', 'C0879590', 'C1334291', 'C1333539', 'C1367721', 'C1334479',
  'C0919466', 'C1335284', 'C1418901', 'C0812215', 'C1970159', 'C1366644',
  'C1333572', 'C1333669', 'C1419044', 'C1420314', 'C1419030', 'C1335212',
  'C1334140', 'C0812265', 'C1332790', 'C1333535', 'C1425013', 'C1334522',
  'C0919507', 'C1414607', 'C1412166', 'C1418926', 'C0080298', 'C1420385',
  'C1414603', 'C1418576', 'C1334117', 'C1336620', 'C1414593', 'C1414592',
  'C1336929', 'C1538740', 'C1335202', 'C1334906', 'C1422166', 'C1332791',
  'C1335201', 'C1334480', 'C1417230', 'C1334863', 'C1425540', 'C1367472',
  'C1366894', 'C1425048', 'C0242957', 'C1335823', 'C1416498', 'C1421267',
  'C1419280', 'C1332453', 'C1333707', 'C1414655', 'C1333928', 'C1418903',
  'C1333578', 'C1366765', 'C1334292', 'C1857782', 'C1333362', 'C1419282',
  'C1417668', 'C1419037', 'C1419041', 'C1333543', 'C1367476', 'C0812290',
  'C1334116'
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
//     --nodes         c4_mapk_nodes.csv \
//     --edges         c4_mapk_edges_all.csv \
//     --mapping       config/bifo_ddkg_mapping.yaml \
//     --seed-nodes    data/benchmark/c4_mapk_seed_nodes.txt \
//     --heldout-nodes data/benchmark/c4_mapk_heldout_nodes.txt \
//     --out-json      results/c4_mapk/results.json
//
//   python score_pathways.py \
//     --nodes             c4_mapk_nodes.csv \
//     --edges-raw         c4_mapk_edges_all.csv \
//     --edges-conditioned results/c4_mapk/results_kept_edges.csv \
//     --scores-cond       results/c4_mapk/results_scores_cond.npy \
//     --scores-raw        results/c4_mapk/results_scores_raw.npy \
//     --node-index        results/c4_mapk/results_node_index.json \
//     --seed-nodes        data/benchmark/c4_mapk_seed_nodes.txt \
//     --chd-pathways      data/benchmark/c4_mapk_pathway_reference.txt \
//     --out-csv           results/c4_mapk/pathway_scores.csv \
//     --out-json          results/c4_mapk/pathway_metrics.json
// =============================================================================
