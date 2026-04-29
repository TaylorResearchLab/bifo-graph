// BIFO CHD Export -- QUERY 5: Pathway membership edges (global)
// Save output as: kf_chd_pathway_membership_edges.csv
// Columns: source, target, predicate
//
// Exports gene->pathway membership edges only (G->PW direction).
// Signal flows FROM genes TO pathways in the BIFO PPR operator.
// Forward direction (PW->G) is nonpropagating_context in the BIFO YAML
// and would be filtered out by conditioning. Excluded at cypher level
// for efficiency.
// Cohort-independent: this query ignores the seed list intentionally.
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
