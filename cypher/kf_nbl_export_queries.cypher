// =============================================================================
// BIFO Kids First Export — Neuroblastoma (KF_NBL)
// Seeds: 88 genes (AutoGVP P/LP, MAF<=0.0001 gnomAD, n_carriers>=3)
// Top genes by carrier count: MAP2K7(12), FLNA(11), ATP1A3(9), GFM1(8), TAF4(8)
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
    "MAP2K7 gene", "FLNA gene", "ATP1A3 gene", "GFM1 gene", "TAF4 gene",
    "ABCB11 gene", "WNK1 gene", "LRP1 gene", "SMPD1 gene", "SLC7A9 gene",
    "GALT gene", "SEC63 gene", "PC gene", "NSD1 gene", "PRMT1 gene",
    "ENO3 gene", "ATRX gene", "KIAA1109 gene", "VWF gene", "SKIV2L gene",
    "GSS gene", "BARD1 gene", "KCNT2 gene", "VCAN gene", "SHMT2 gene",
    "PIEZO1 gene", "PCK1 gene", "FH gene", "FANCA gene", "PDE6B gene",
    "CHRNE gene", "PAH gene", "SLC12A3 gene", "BLM gene", "ZSWIM6 gene",
    "TYR gene", "ANO3 gene", "DYNC2H1 gene", "POU6F2 gene", "ARMC5 gene",
    "STAG3 gene", "CDC45 gene", "ABCB6 gene", "PGD gene", "NCAPD3 gene",
    "FSIP2 gene", "CLPP gene", "CFTR gene", "SEPSECS gene", "ATM gene",
    "ULK1 gene", "HLA-DRB1 gene", "MAF gene", "DNAH5 gene", "SLC45A2 gene",
    "ZFX gene", "COL7A1 gene", "DHCR24 gene", "DNAH9 gene", "GBE1 gene",
    "SPG7 gene", "PLXND1 gene", "HPS3 gene", "AHCY gene", "TTC21B gene",
    "RAD54B gene", "DUOXA2 gene", "ERCC2 gene", "GRM6 gene", "FREM1 gene",
    "HK1 gene", "MED12 gene", "LONP1 gene", "G6PC1 gene", "NIN gene",
    "ADGRV1 gene", "CHEK2 gene", "NUP205 gene", "OGDH gene", "WFS1 gene",
    "IRF2BPL gene", "HNF1B gene", "ASPM gene", "VPS13B gene", "HSPA9 gene",
    "NCOR1 gene", "MAT1A gene", "ABCA13 gene"
] AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
RETURN replace(t.name, " gene", "") AS gene_symbol,
       c.CUI                        AS cui,
       t.name                       AS pref_term
ORDER BY gene_symbol;


// ---------------------------------------------------------------------------
// QUERY 2 — Seed concept nodes → save as: kf_nbl_seed_nodes.csv
// Columns: node_id, label, name, sab
// ---------------------------------------------------------------------------
WITH [
    "MAP2K7 gene", "FLNA gene", "ATP1A3 gene", "GFM1 gene", "TAF4 gene",
    "ABCB11 gene", "WNK1 gene", "LRP1 gene", "SMPD1 gene", "SLC7A9 gene",
    "GALT gene", "SEC63 gene", "PC gene", "NSD1 gene", "PRMT1 gene",
    "ENO3 gene", "ATRX gene", "KIAA1109 gene", "VWF gene", "SKIV2L gene",
    "GSS gene", "BARD1 gene", "KCNT2 gene", "VCAN gene", "SHMT2 gene",
    "PIEZO1 gene", "PCK1 gene", "FH gene", "FANCA gene", "PDE6B gene",
    "CHRNE gene", "PAH gene", "SLC12A3 gene", "BLM gene", "ZSWIM6 gene",
    "TYR gene", "ANO3 gene", "DYNC2H1 gene", "POU6F2 gene", "ARMC5 gene",
    "STAG3 gene", "CDC45 gene", "ABCB6 gene", "PGD gene", "NCAPD3 gene",
    "FSIP2 gene", "CLPP gene", "CFTR gene", "SEPSECS gene", "ATM gene",
    "ULK1 gene", "HLA-DRB1 gene", "MAF gene", "DNAH5 gene", "SLC45A2 gene",
    "ZFX gene", "COL7A1 gene", "DHCR24 gene", "DNAH9 gene", "GBE1 gene",
    "SPG7 gene", "PLXND1 gene", "HPS3 gene", "AHCY gene", "TTC21B gene",
    "RAD54B gene", "DUOXA2 gene", "ERCC2 gene", "GRM6 gene", "FREM1 gene",
    "HK1 gene", "MED12 gene", "LONP1 gene", "G6PC1 gene", "NIN gene",
    "ADGRV1 gene", "CHEK2 gene", "NUP205 gene", "OGDH gene", "WFS1 gene",
    "IRF2BPL gene", "HNF1B gene", "ASPM gene", "VPS13B gene", "HSPA9 gene",
    "NCOR1 gene", "MAT1A gene", "ABCA13 gene"
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
// QUERY 3 — 1-hop edges → save as: kf_nbl_edges_raw.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
WITH [
    "MAP2K7 gene", "FLNA gene", "ATP1A3 gene", "GFM1 gene", "TAF4 gene",
    "ABCB11 gene", "WNK1 gene", "LRP1 gene", "SMPD1 gene", "SLC7A9 gene",
    "GALT gene", "SEC63 gene", "PC gene", "NSD1 gene", "PRMT1 gene",
    "ENO3 gene", "ATRX gene", "KIAA1109 gene", "VWF gene", "SKIV2L gene",
    "GSS gene", "BARD1 gene", "KCNT2 gene", "VCAN gene", "SHMT2 gene",
    "PIEZO1 gene", "PCK1 gene", "FH gene", "FANCA gene", "PDE6B gene",
    "CHRNE gene", "PAH gene", "SLC12A3 gene", "BLM gene", "ZSWIM6 gene",
    "TYR gene", "ANO3 gene", "DYNC2H1 gene", "POU6F2 gene", "ARMC5 gene",
    "STAG3 gene", "CDC45 gene", "ABCB6 gene", "PGD gene", "NCAPD3 gene",
    "FSIP2 gene", "CLPP gene", "CFTR gene", "SEPSECS gene", "ATM gene",
    "ULK1 gene", "HLA-DRB1 gene", "MAF gene", "DNAH5 gene", "SLC45A2 gene",
    "ZFX gene", "COL7A1 gene", "DHCR24 gene", "DNAH9 gene", "GBE1 gene",
    "SPG7 gene", "PLXND1 gene", "HPS3 gene", "AHCY gene", "TTC21B gene",
    "RAD54B gene", "DUOXA2 gene", "ERCC2 gene", "GRM6 gene", "FREM1 gene",
    "HK1 gene", "MED12 gene", "LONP1 gene", "G6PC1 gene", "NIN gene",
    "ADGRV1 gene", "CHEK2 gene", "NUP205 gene", "OGDH gene", "WFS1 gene",
    "IRF2BPL gene", "HNF1B gene", "ASPM gene", "VPS13B gene", "HSPA9 gene",
    "NCOR1 gene", "MAT1A gene", "ABCA13 gene"
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
// QUERY 4 — All nodes (seeds + neighbors) → save as: kf_nbl_nodes.csv
// Columns: node_id, label, name, sab
// Uses CASE to prefer HGNC > NCBIGENE > MSIGDB > GO > MONDO > OMIM
// This is critical — bifo_conditioning.py uses sab for entity resolution
// ---------------------------------------------------------------------------
WITH [
    "MAP2K7 gene", "FLNA gene", "ATP1A3 gene", "GFM1 gene", "TAF4 gene",
    "ABCB11 gene", "WNK1 gene", "LRP1 gene", "SMPD1 gene", "SLC7A9 gene",
    "GALT gene", "SEC63 gene", "PC gene", "NSD1 gene", "PRMT1 gene",
    "ENO3 gene", "ATRX gene", "KIAA1109 gene", "VWF gene", "SKIV2L gene",
    "GSS gene", "BARD1 gene", "KCNT2 gene", "VCAN gene", "SHMT2 gene",
    "PIEZO1 gene", "PCK1 gene", "FH gene", "FANCA gene", "PDE6B gene",
    "CHRNE gene", "PAH gene", "SLC12A3 gene", "BLM gene", "ZSWIM6 gene",
    "TYR gene", "ANO3 gene", "DYNC2H1 gene", "POU6F2 gene", "ARMC5 gene",
    "STAG3 gene", "CDC45 gene", "ABCB6 gene", "PGD gene", "NCAPD3 gene",
    "FSIP2 gene", "CLPP gene", "CFTR gene", "SEPSECS gene", "ATM gene",
    "ULK1 gene", "HLA-DRB1 gene", "MAF gene", "DNAH5 gene", "SLC45A2 gene",
    "ZFX gene", "COL7A1 gene", "DHCR24 gene", "DNAH9 gene", "GBE1 gene",
    "SPG7 gene", "PLXND1 gene", "HPS3 gene", "AHCY gene", "TTC21B gene",
    "RAD54B gene", "DUOXA2 gene", "ERCC2 gene", "GRM6 gene", "FREM1 gene",
    "HK1 gene", "MED12 gene", "LONP1 gene", "G6PC1 gene", "NIN gene",
    "ADGRV1 gene", "CHEK2 gene", "NUP205 gene", "OGDH gene", "WFS1 gene",
    "IRF2BPL gene", "HNF1B gene", "ASPM gene", "VPS13B gene", "HSPA9 gene",
    "NCOR1 gene", "MAT1A gene", "ABCA13 gene"
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
// QUERY 5 — Pathway membership edges → save as: kf_nbl_pathway_membership_edges.csv
// Columns: source, target, predicate
// ---------------------------------------------------------------------------
WITH [
    "MAP2K7 gene", "FLNA gene", "ATP1A3 gene", "GFM1 gene", "TAF4 gene",
    "ABCB11 gene", "WNK1 gene", "LRP1 gene", "SMPD1 gene", "SLC7A9 gene",
    "GALT gene", "SEC63 gene", "PC gene", "NSD1 gene", "PRMT1 gene",
    "ENO3 gene", "ATRX gene", "KIAA1109 gene", "VWF gene", "SKIV2L gene",
    "GSS gene", "BARD1 gene", "KCNT2 gene", "VCAN gene", "SHMT2 gene",
    "PIEZO1 gene", "PCK1 gene", "FH gene", "FANCA gene", "PDE6B gene",
    "CHRNE gene", "PAH gene", "SLC12A3 gene", "BLM gene", "ZSWIM6 gene",
    "TYR gene", "ANO3 gene", "DYNC2H1 gene", "POU6F2 gene", "ARMC5 gene",
    "STAG3 gene", "CDC45 gene", "ABCB6 gene", "PGD gene", "NCAPD3 gene",
    "FSIP2 gene", "CLPP gene", "CFTR gene", "SEPSECS gene", "ATM gene",
    "ULK1 gene", "HLA-DRB1 gene", "MAF gene", "DNAH5 gene", "SLC45A2 gene",
    "ZFX gene", "COL7A1 gene", "DHCR24 gene", "DNAH9 gene", "GBE1 gene",
    "SPG7 gene", "PLXND1 gene", "HPS3 gene", "AHCY gene", "TTC21B gene",
    "RAD54B gene", "DUOXA2 gene", "ERCC2 gene", "GRM6 gene", "FREM1 gene",
    "HK1 gene", "MED12 gene", "LONP1 gene", "G6PC1 gene", "NIN gene",
    "ADGRV1 gene", "CHEK2 gene", "NUP205 gene", "OGDH gene", "WFS1 gene",
    "IRF2BPL gene", "HNF1B gene", "ASPM gene", "VPS13B gene", "HSPA9 gene",
    "NCOR1 gene", "MAT1A gene", "ABCA13 gene"
] AS gene_names

MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
WITH collect(DISTINCT c) AS seeds

UNWIND seeds AS seed
OPTIONAL MATCH (seed)-[r1]-(neighbor:Concept)
WHERE type(r1) <> 'CODE'
  AND type(r1) <> 'STY'
  AND type(r1) <> 'ISA'
  AND type(r1) <> 'PREF_TERM'
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
// AFTER EXPORT — merge with pandas (not cat):
//
//   python3 -c "
//   import pandas as pd
//   pd.concat([pd.read_csv('kf_nbl_nodes.csv'), pd.read_csv('kf_nbl_ncc_pathway_nodes.csv')],
//             ignore_index=True).to_csv('kf_nbl_nodes_extended.csv', index=False)
//   pd.concat([pd.read_csv('kf_nbl_edges_raw.csv'), pd.read_csv('kf_nbl_pathway_membership_edges.csv')],
//             ignore_index=True).to_csv('kf_nbl_edges_merged.csv', index=False)
//   pd.concat([pd.read_csv('kf_nbl_edges_merged.csv'), pd.read_csv('kf_nbl_ncc_membership_edges.csv')],
//             ignore_index=True).to_csv('kf_nbl_edges_all.csv', index=False)
//   "
// =============================================================================
