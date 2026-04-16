// =============================================================================
// BIFO CHD Curated Benchmark — Neo4j Export Queries
// (originally: minimal_test_run_neo4j_export.cypher)
// =============================================================================
// Exports a 1-hop concept neighborhood around 15 CHD genes (10 seeds + 5 heldout).
// Run each query separately in Neo4j Browser or cypher-shell.
// Output: nodes.csv and edges_raw.csv
//
// IMPORTANT: Exported CSV files must have these columns:
//   nodes.csv     — node_id, label, name, sab
//   edges_raw.csv — source, target, predicate
//
// The 'sab' column on nodes.csv is required by score_pathways.py.
// The 'label' column must contain 'Concept' for concept nodes.
// =============================================================================

// --- Step 1: Verify seed CUIs (already done — results below for reference) ---
// NOTE: HGNC preferred terms in this DDKG build are stored as
// "SYMBOL gene" (e.g. "GATA4 gene"), not bare symbols.
//
// Verified CUIs (from live DDKG query):
// FLT4    -> C1333569  HGNC:3767   (heldout)
// GATA4   -> C1414995  HGNC:4173   (seed)
// GATA6   -> C1414996  HGNC:4174   (seed)
// HAND1   -> C1415466  HGNC:4807   (seed)
// HAND2   -> C1415467  HGNC:4808   (seed)
// JAG1    -> C1416525  HGNC:6188   (heldout)
// MYH6    -> C1417541  HGNC:7576   (seed)
// MYH7    -> C1417542  HGNC:7577   (heldout)
// NKX2-5  -> C1413784  HGNC:2488   (seed)
// NOTCH1  -> C1334889  HGNC:7881   (seed)
// NOTCH2  -> C1417767  HGNC:7882   (seed)
// PTPN11  -> C1335280  HGNC:9644   (heldout)
// TBX1    -> C1420603  HGNC:11592  (seed)
// TBX5    -> C1420615  HGNC:11604  (seed)
// ZFPM2   -> C1424490  HGNC:16700  (heldout)
//
// To re-run Step 1:
// MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
// MATCH (c)-[:CODE]->(code:Code)
// WHERE code.SAB = 'HGNC'
//   AND t.name IN [
//     'GATA4 gene', 'NKX2-5 gene', 'TBX5 gene', 'NOTCH1 gene', 'NOTCH2 gene',
//     'HAND1 gene', 'HAND2 gene', 'MYH6 gene', 'GATA6 gene', 'TBX1 gene',
//     'ZFPM2 gene', 'MYH7 gene', 'PTPN11 gene', 'JAG1 gene', 'FLT4 gene'
//   ]
// RETURN t.name AS gene, c.CUI AS cui, code.CodeID AS hgnc_id
// ORDER BY gene;

// =============================================================================
// --- Step 2: Export 1-hop neighborhood nodes ---
// Save result as nodes.csv
// Expected scale: ~500–5,000 rows
// =============================================================================

MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1333569', 'C1414995', 'C1414996', 'C1415466', 'C1415467',
  'C1416525', 'C1417541', 'C1417542', 'C1413784', 'C1334889',
  'C1417767', 'C1335280', 'C1420603', 'C1420615', 'C1424490'
]
MATCH (seed)-[r1]-(hop1:Concept)
WHERE type(r1) <> 'CODE' AND type(r1) <> 'PREF_TERM'

WITH collect(DISTINCT seed) + collect(DISTINCT hop1) AS all_concepts

UNWIND all_concepts AS c
WITH DISTINCT c

// Get one representative SAB per Concept (prefer HGNC > NCBIGENE > MSIGDB > GO > MONDO)
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

// Use :Term label (confirmed from diagnostic query — PREF_TERM -> Term nodes)
OPTIONAL MATCH (c)-[:PREF_TERM]->(pt:Term)

RETURN
  c.CUI                              AS node_id,
  'Concept'                          AS label,
  coalesce(pt.name, c.CUI)           AS name,
  coalesce(best_code.SAB, 'UNKNOWN') AS sab;

// =============================================================================
// --- Step 3: Export edges (seed ↔ hop1 only) ---
// Save result as edges_raw.csv
// Expected scale: ~5,000–50,000 rows
//
// NOTE: This query collects only seed↔hop1 edges in a single traversal pass.
// It deliberately excludes hop1↔hop1 edges — those would pull in massive
// cross-talk between unrelated neighbors (coexpression, clinical coding,
// drug relationships) not driven by CHD signal.
// For PPR purposes this is correct: signal propagates from seeds into hop1
// nodes via these edges; lateral hop1↔hop1 noise is excluded.
// =============================================================================

MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1333569', 'C1414995', 'C1414996', 'C1415466', 'C1415467',
  'C1416525', 'C1417541', 'C1417542', 'C1413784', 'C1334889',
  'C1417767', 'C1335280', 'C1420603', 'C1420615', 'C1424490'
]
MATCH (seed)-[r]-(hop1:Concept)
WHERE type(r) <> 'CODE' AND type(r) <> 'PREF_TERM'
RETURN
  startNode(r).CUI AS source,
  endNode(r).CUI   AS target,
  type(r)          AS predicate;

// =============================================================================
// --- Step 4: Export pathway membership edges (targeted hop1↔hop1) ---
// Save result as pathway_membership_edges.csv
// Then append to edges_raw.csv before running the pipeline:
//   tail -n +2 pathway_membership_edges.csv >> edges_raw.csv
//
// WHY: The Step 3 export only includes seed↔hop1 edges.
// Pathway membership edges (MSIGDB→HGNC gene) connect two hop1 nodes —
// they are absent from edges_raw.csv. This step adds only the specific
// membership predicates needed for member-gene aggregation in score_pathways.py.
// Using a targeted predicate list avoids re-opening the full hop1↔hop1 flood.
// =============================================================================

MATCH (seed:Concept)
WHERE seed.CUI IN [
  'C1333569', 'C1414995', 'C1414996', 'C1415466', 'C1415467',
  'C1416525', 'C1417541', 'C1417542', 'C1413784', 'C1334889',
  'C1417767', 'C1335280', 'C1420603', 'C1420615', 'C1424490'
]
MATCH (seed)-[r1]-(hop1:Concept)
WHERE type(r1) <> 'CODE' AND type(r1) <> 'PREF_TERM'
WITH collect(DISTINCT hop1.CUI) AS hop1_ids

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

RETURN
  startNode(r).CUI AS source,
  endNode(r).CUI   AS target,
  type(r)          AS predicate;
