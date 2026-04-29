// BIFO CHD Export -- QUERY 6: Pathway member nodes (full MSIGDB membership)
// Save output as: kf_chd_pathway_member_nodes.csv
// Columns: node_id, label, name, sab
//
// Exports node metadata for all genes that are MSIGDB pathway members.
// These genes may not appear in the 1-hop seed neighborhood but are required
// for entity resolution in bifo_conditioning.py when their membership edges
// are present in the pathway_membership_edges file.
// Cohort-independent: this query ignores the seed list intentionally.
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
