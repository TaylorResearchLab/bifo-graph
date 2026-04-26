// =============================================================================
// BIFO Kids First Export — Congenital Heart Defects (KF_CHD)
// Seeds: 56 genes (AutoGVP P/LP, MAF<=0.0001 gnomAD, n_carriers>=3)
// Top genes by carrier count: ROS1(5), SLC12A3(5), VPS11(5), DNAH5(5), HPS3(5)
//
// Schema: (c:Concept)-[:PREF_TERM]->(t:Term {name:"SYMBOL gene"})
// Output columns match bifo_conditioning.py exactly:
//   nodes    -> node_id, label, name, sab
//   edges    -> source, target, predicate
// =============================================================================

// ---------------------------------------------------------------------------
// QUERY 1 — VERIFY gene coverage (run in Browser, no CSV needed)
// ---------------------------------------------------------------------------
WITH [
    "ROS1 gene", "SLC12A3 gene", "VPS11 gene", "DNAH5 gene", "HPS3 gene",
    "BRCA2 gene", "CC2D2A gene", "SPG11 gene", "GPI gene", "DUOX2 gene",
    "ACADVL gene", "ABCA13 gene", "ASPM gene", "CEP290 gene", "NOTCH2 gene",
    "FANCA gene", "AGL gene", "CACNA1S gene", "KIAA0586 gene", "FANCI gene",
    "GBE1 gene", "DZIP1L gene", "LRBA gene", "OTOF gene", "TTI1 gene",
    "SYNE1 gene", "ABCA7 gene", "CLCN3 gene", "AAAS gene", "PDZRN3 gene",
    "PKD1L1 gene", "TRNT1 gene", "PDHA2 gene", "ACE gene", "SLC26A4 gene",
    "LDLR gene", "LAMA5 gene", "MEI1 gene", "SLCO1B1 gene", "CUBN gene",
    "NEB gene", "THRB gene", "TMC1 gene", "DHCR7 gene", "CD36 gene",
    "DMBX1 gene", "PTPRQ gene", "CERKL gene", "SLC17A9 gene", "PIK3R5 gene",
    "ARSA gene", "MAP3K6 gene", "SLC24A1 gene", "ADCY10 gene", "ATP7B gene",
    "DRD4 gene"
] AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
RETURN replace(t.name, " gene", "") AS gene_symbol,
       c.CUI                        AS cui,
       t.name                       AS pref_term
ORDER BY gene_symbol;


// ---------------------------------------------------------------------------
// QUERY 2 — Seed concept nodes → save as: kf_chd_seed_nodes.csv
// Columns: node_id, label, name, sab
// ---------------------------------------------------------------------------
WITH [
    "ROS1 gene", "SLC12A3 gene", "VPS11 gene", "DNAH5 gene", "HPS3 gene",
    "BRCA2 gene", "CC2D2A gene", "SPG11 gene", "GPI gene", "DUOX2 gene",
    "ACADVL gene", "ABCA13 gene", "ASPM gene", "CEP290 gene", "NOTCH2 gene",
    "FANCA gene", "AGL gene", "CACNA1S gene", "KIAA0586 gene", "FANCI gene",
    "GBE1 gene", "DZIP1L gene", "LRBA gene", "OTOF gene", "TTI1 gene",
    "SYNE1 gene", "ABCA7 gene", "CLCN3 gene", "AAAS gene", "PDZRN3 gene",
    "PKD1L1 gene", "TRNT1 gene", "PDHA2 gene", "ACE gene", "SLC26A4 gene",
    "LDLR gene", "LAMA5 gene", "MEI1 gene", "SLCO1B1 gene", "CUBN gene",
    "NEB gene", "THRB gene", "TMC1 gene", "DHCR7 gene", "CD36 gene",
    "DMBX1 gene", "PTPRQ gene", "CERKL gene", "SLC17A9 gene", "PIK3R5 gene",
    "ARSA gene", "MAP3K6 gene", "SLC24A1 gene", "ADCY10 gene", "ATP7B gene",
    "DRD4 gene"
] AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
OPTIONAL MATCH (c)-[:CODE]->(code:Code)
WITH c, t, code
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
WITH c, t, collect(code)[0] AS best_code
RETURN DISTINCT
    c.CUI                              AS node_id,
    'Concept'                          AS label,
    coalesce(t.name, c.CUI)            AS name,
    coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;


// ---------------------------------------------------------------------------
// QUERY 3 — 1-hop edges → save as: kf_chd_edges_raw.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
WITH [
    "ROS1 gene", "SLC12A3 gene", "VPS11 gene", "DNAH5 gene", "HPS3 gene",
    "BRCA2 gene", "CC2D2A gene", "SPG11 gene", "GPI gene", "DUOX2 gene",
    "ACADVL gene", "ABCA13 gene", "ASPM gene", "CEP290 gene", "NOTCH2 gene",
    "FANCA gene", "AGL gene", "CACNA1S gene", "KIAA0586 gene", "FANCI gene",
    "GBE1 gene", "DZIP1L gene", "LRBA gene", "OTOF gene", "TTI1 gene",
    "SYNE1 gene", "ABCA7 gene", "CLCN3 gene", "AAAS gene", "PDZRN3 gene",
    "PKD1L1 gene", "TRNT1 gene", "PDHA2 gene", "ACE gene", "SLC26A4 gene",
    "LDLR gene", "LAMA5 gene", "MEI1 gene", "SLCO1B1 gene", "CUBN gene",
    "NEB gene", "THRB gene", "TMC1 gene", "DHCR7 gene", "CD36 gene",
    "DMBX1 gene", "PTPRQ gene", "CERKL gene", "SLC17A9 gene", "PIK3R5 gene",
    "ARSA gene", "MAP3K6 gene", "SLC24A1 gene", "ADCY10 gene", "ATP7B gene",
    "DRD4 gene"
] AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
WITH collect(DISTINCT c) AS seed_concepts

UNWIND seed_concepts AS seed
MATCH (seed)-[r]-(neighbor:Concept)
WHERE type(r) <> 'CODE'
  AND type(r) <> 'STY'
  AND type(r) <> 'ISA'
  AND type(r) <> 'PREF_TERM'
RETURN DISTINCT
    seed.CUI     AS source,
    neighbor.CUI AS target,
    type(r)      AS predicate
ORDER BY source, predicate, target;


// ---------------------------------------------------------------------------
// QUERY 4 — All nodes (seeds + neighbors) → save as: kf_chd_nodes.csv
// Columns: node_id, label, name, sab
// Uses CASE to prefer HGNC > NCBIGENE > MSIGDB > GO > MONDO > OMIM
// This is critical — bifo_conditioning.py uses sab for entity resolution
// ---------------------------------------------------------------------------
WITH [
    "ROS1 gene", "SLC12A3 gene", "VPS11 gene", "DNAH5 gene", "HPS3 gene",
    "BRCA2 gene", "CC2D2A gene", "SPG11 gene", "GPI gene", "DUOX2 gene",
    "ACADVL gene", "ABCA13 gene", "ASPM gene", "CEP290 gene", "NOTCH2 gene",
    "FANCA gene", "AGL gene", "CACNA1S gene", "KIAA0586 gene", "FANCI gene",
    "GBE1 gene", "DZIP1L gene", "LRBA gene", "OTOF gene", "TTI1 gene",
    "SYNE1 gene", "ABCA7 gene", "CLCN3 gene", "AAAS gene", "PDZRN3 gene",
    "PKD1L1 gene", "TRNT1 gene", "PDHA2 gene", "ACE gene", "SLC26A4 gene",
    "LDLR gene", "LAMA5 gene", "MEI1 gene", "SLCO1B1 gene", "CUBN gene",
    "NEB gene", "THRB gene", "TMC1 gene", "DHCR7 gene", "CD36 gene",
    "DMBX1 gene", "PTPRQ gene", "CERKL gene", "SLC17A9 gene", "PIK3R5 gene",
    "ARSA gene", "MAP3K6 gene", "SLC24A1 gene", "ADCY10 gene", "ATP7B gene",
    "DRD4 gene"
] AS gene_names

MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
WITH collect(DISTINCT c) AS seed_concepts

UNWIND seed_concepts AS seed
OPTIONAL MATCH (seed)-[r]-(neighbor:Concept)
WHERE type(r) <> 'CODE'
  AND type(r) <> 'STY'
  AND type(r) <> 'ISA'
  AND type(r) <> 'PREF_TERM'
WITH seed_concepts, collect(DISTINCT neighbor) AS neighbors
WITH seed_concepts + neighbors AS all_concepts

UNWIND all_concepts AS c
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
RETURN DISTINCT
    c.CUI                              AS node_id,
    'Concept'                          AS label,
    coalesce(pt.name, c.CUI)           AS name,
    coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// QUERY 5 — MSIGDB gene→pathway membership edges → save as: kf_chd_pathway_membership_edges.csv
// Columns: source, target, predicate
//
// Exports gene→pathway membership edges only (unidirectional).
// Signal flows FROM genes TO pathways in the BIFO PPR operator.
// Pathway→gene direction (pathway_associated_with_gene) is excluded to
// prevent feedback amplification loops in the PPR operator.
// Includes both MSIGDB C2.CP (inverse_pathway_associated_with_gene)
// and Hallmark (inverse_has_signature_gene) membership predicates.
// ---------------------------------------------------------------------------
MATCH (pw:Concept)-[:CODE]->(code:Code)
WHERE code.SAB = 'MSIGDB'
WITH pw
MATCH (gene:Concept)-[:inverse_pathway_associated_with_gene]->(pw)
RETURN DISTINCT
    gene.CUI AS source,
    pw.CUI   AS target,
    'inverse_pathway_associated_with_gene' AS predicate
UNION
MATCH (pw:Concept)-[:CODE]->(code:Code)
WHERE code.SAB = 'MSIGDB'
WITH pw
MATCH (gene:Concept)-[:inverse_has_signature_gene]->(pw)
RETURN DISTINCT
    gene.CUI AS source,
    pw.CUI   AS target,
    'inverse_has_signature_gene' AS predicate
ORDER BY source, target;


// ---------------------------------------------------------------------------
// QUERY 6 — Gene nodes for full MSIGDB membership → save as: kf_chd_pathway_member_nodes.csv
// Columns: node_id, label, name, sab
//
// Exports node metadata for all genes that are MSIGDB pathway members.
// These genes may not appear in the 1-hop seed neighborhood but are required
// for entity resolution in bifo_conditioning.py when their membership edges
// are present in the pathway_membership_edges file.
// SAB priority matches Query 4: HGNC > NCBIGENE > MSIGDB > GO > MONDO > OMIM
// ---------------------------------------------------------------------------
MATCH (pw:Concept)-[:CODE]->(code:Code)
WHERE code.SAB = 'MSIGDB'
WITH pw
MATCH (pw)-[:pathway_associated_with_gene]->(gene:Concept)
WITH collect(DISTINCT gene) AS all_genes
UNWIND all_genes AS c
OPTIONAL MATCH (c)-[:CODE]->(code:Code)
WITH c, code
ORDER BY
  CASE code.SAB
    WHEN 'HGNC'     THEN 0
    WHEN 'NCBIGENE' THEN 1
    WHEN 'MSIGDB'   THEN 2
    WHEN 'GO'       THEN 3
    WHEN 'MONDO'    THEN 4
    WHEN 'OMIM'     THEN 5
    ELSE 99
  END
WITH c, collect(code)[0] AS best_code
OPTIONAL MATCH (c)-[:PREF_TERM]->(t:Term)
RETURN DISTINCT
    c.CUI                              AS node_id,
    'Concept'                          AS label,
    coalesce(t.name, c.CUI)            AS name,
    coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;


// =============================================================================
// =============================================================================
// AFTER EXPORT — merge with pandas (not cat):
//
//   python3 -c "
//   import pandas as pd
//   pd.concat([pd.read_csv('kf_chd_nodes.csv'), pd.read_csv('kf_chd_ncc_pathway_nodes.csv')],
//             ignore_index=True).to_csv('kf_chd_nodes_extended.csv', index=False)
//   pd.concat([pd.read_csv('kf_chd_edges_raw.csv'), pd.read_csv('kf_chd_pathway_membership_edges.csv')],
//             ignore_index=True).to_csv('kf_chd_edges_merged.csv', index=False)
//   pd.concat([pd.read_csv('kf_chd_edges_merged.csv'), pd.read_csv('kf_chd_ncc_membership_edges.csv')],
//             ignore_index=True).to_csv('kf_chd_edges_all.csv', index=False)
//   "
// =============================================================================
