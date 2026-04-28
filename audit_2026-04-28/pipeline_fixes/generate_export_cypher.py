#!/usr/bin/env python3
"""
generate_export_cypher.py
Generates cypher-shell query files (query2-6) for BIFO Neo4j export.

Reads a seed gene list (HGNC symbols) and produces five split .cypher files
ready to be run by run_kf_{cohort}_export.sh via cypher-shell.

Usage:
    python generate_export_cypher.py \
        --seeds   kf_chd_seeds_maf001.txt \
        --cohort  chd \
        --out-dir .

Output files (written to --out-dir):
    kf_{cohort}_query2.cypher   -- seed concept nodes      -> kf_{cohort}_seed_nodes.csv
    kf_{cohort}_query3.cypher   -- 1-hop edges             -> kf_{cohort}_edges_raw.csv
    kf_{cohort}_query4.cypher   -- all nodes               -> kf_{cohort}_nodes.csv
    kf_{cohort}_query5.cypher   -- pathway membership      -> kf_{cohort}_pathway_membership_edges.csv
    kf_{cohort}_query6.cypher   -- pathway member nodes    -> kf_{cohort}_pathway_member_nodes.csv

CHANGELOG (audit_2026-04-28):
    Added make_query6() (ported from examples/strict_filter/) to produce the
    pathway_member_nodes.csv file expected by clean_files.sh and merge_files.sh.
    The production cypher (committed in 5d01d32) had only queries 2-5; that gap
    meant the documented pipeline could not run end-to-end without manual
    intervention. Q6 exports node metadata for all genes that are MSIGDB pathway
    members, required for entity resolution in bifo_conditioning.py when
    membership edges reference genes outside the 1-hop seed neighborhood.
"""

import argparse
from pathlib import Path


def load_genes(seeds_file):
    genes = []
    for line in Path(seeds_file).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        gene = line.split('\t')[0].strip()
        if gene:
            genes.append(gene)
    return genes


def format_gene_list(genes, indent=4):
    """Format gene list as Cypher string array, 5 per line."""
    pad = ' ' * indent
    quoted = [f'"{g} gene"' for g in genes]
    lines = []
    for i in range(0, len(quoted), 5):
        lines.append(pad + ', '.join(quoted[i:i+5]))
    return '[\n' + ',\n'.join(lines) + '\n]'


SAB_CASE = """  CASE code.SAB
    WHEN 'HGNC'     THEN 0
    WHEN 'NCBIGENE' THEN 1
    WHEN 'MSIGDB'   THEN 2
    WHEN 'REACTOME' THEN 3
    WHEN 'GO'       THEN 4
    WHEN 'MONDO'    THEN 5
    WHEN 'OMIM'     THEN 6
    ELSE 99
  END ASC"""


def make_query2(cohort, gene_list_str, n_genes):
    return f"""// BIFO {cohort.upper()} Export -- QUERY 2: Seed concept nodes
// {n_genes} seed genes (AutoGVP P/LP)
// Save output as: kf_{cohort}_seed_nodes.csv
// Columns: node_id, label, name, sab
WITH {gene_list_str} AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
OPTIONAL MATCH (c)-[:CODE]->(code:Code)
WITH c, t, code
ORDER BY
{SAB_CASE}
WITH c, t, collect(code)[0] AS best_code
RETURN DISTINCT
    c.CUI                              AS node_id,
    'Concept'                          AS label,
    coalesce(t.name, c.CUI)            AS name,
    coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;
"""


def make_query3(cohort, gene_list_str, n_genes):
    return f"""// BIFO {cohort.upper()} Export -- QUERY 3: 1-hop edges
// Save output as: kf_{cohort}_edges_raw.csv
// Columns: source, target, predicate
WITH {gene_list_str} AS gene_names
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
"""


def make_query4(cohort, gene_list_str, n_genes):
    return f"""// BIFO {cohort.upper()} Export -- QUERY 4: All nodes (seeds + 1-hop neighbors)
// Save output as: kf_{cohort}_nodes.csv
// Columns: node_id, label, name, sab
WITH {gene_list_str} AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
WITH collect(DISTINCT c) AS seed_concepts

UNWIND seed_concepts AS seed
MATCH (seed)-[r]-(hop1:Concept)
WHERE type(r) <> 'CODE'
  AND type(r) <> 'STY'
  AND type(r) <> 'ISA'
  AND type(r) <> 'PREF_TERM'

WITH collect(DISTINCT seed) + collect(DISTINCT hop1) AS all_concepts
UNWIND all_concepts AS c
WITH DISTINCT c

OPTIONAL MATCH (c)-[:CODE]->(code:Code)
WITH c, code
ORDER BY
{SAB_CASE}
WITH c, collect(code)[0] AS best_code
OPTIONAL MATCH (c)-[:PREF_TERM]->(pt:Term)

RETURN
    c.CUI                              AS node_id,
    'Concept'                          AS label,
    coalesce(pt.name, c.CUI)           AS name,
    coalesce(best_code.SAB, 'UNKNOWN') AS sab
ORDER BY node_id;
"""


def make_query5(cohort, gene_list_str, n_genes):
    return f"""// BIFO {cohort.upper()} Export -- QUERY 5: Pathway membership edges
// Save output as: kf_{cohort}_pathway_membership_edges.csv
// Columns: source, target, predicate
WITH {gene_list_str} AS gene_names
MATCH (c:Concept)-[:PREF_TERM]->(t:Term)
WHERE t.name IN gene_names
WITH collect(DISTINCT c) AS seed_concepts

UNWIND seed_concepts AS seed
MATCH (seed)-[r1]-(hop1:Concept)
WHERE type(r1) <> 'CODE'
  AND type(r1) <> 'STY'
  AND type(r1) <> 'ISA'
  AND type(r1) <> 'PREF_TERM'
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
RETURN DISTINCT
    pw.CUI   AS source,
    gene.CUI AS target,
    type(r)  AS predicate
ORDER BY source, predicate, target;
"""


def make_query6(cohort, gene_list_str, n_genes):
    """
    Pathway member nodes: exports node metadata for all genes that are MSIGDB
    pathway members, regardless of whether they appear in the 1-hop seed
    neighborhood. Required for entity resolution in bifo_conditioning.py when
    membership edges reference genes outside the seed neighborhood.

    Ported from examples/strict_filter/kf_chd_export_queries.cypher (audit_2026-04-28).
    """
    return f"""// BIFO {cohort.upper()} Export -- QUERY 6: Pathway member nodes (full MSIGDB membership)
// Save output as: kf_{cohort}_pathway_member_nodes.csv
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
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--seeds',   required=True, help='Seed gene list (.txt)')
    ap.add_argument('--cohort',  required=True,
                    help='Cohort name: chd, nbl, tall, cbtn, rsbd (or other)')
    ap.add_argument('--out-dir', default='.', help='Output directory (default: .)')
    args = ap.parse_args()

    genes = load_genes(args.seeds)
    if not genes:
        raise ValueError(f"No genes found in {args.seeds}")

    gene_list_str = format_gene_list(genes)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    queries = {
        'query2': make_query2(args.cohort, gene_list_str, len(genes)),
        'query3': make_query3(args.cohort, gene_list_str, len(genes)),
        'query4': make_query4(args.cohort, gene_list_str, len(genes)),
        'query5': make_query5(args.cohort, gene_list_str, len(genes)),
        'query6': make_query6(args.cohort, gene_list_str, len(genes)),
    }

    for name, content in queries.items():
        path = out / f'kf_{args.cohort}_{name}.cypher'
        path.write_text(content)
        print(f"Written: {path}")

    print(f"\nGenerated {len(queries)} query files for {args.cohort.upper()}")
    print(f"  {len(genes)} seed genes")
    print(f"  Run: bash run_kf_{args.cohort}_export.sh [user] [pass] [addr]")


if __name__ == '__main__':
    main()
