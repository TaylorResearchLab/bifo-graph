#!/usr/bin/env python3
"""
summarize_results.py — Post-processing script for BIFO pathway analysis results.

Produces two output files per run:

  1. pathway_results_summary.tsv
     Machine-readable table of all scored pathways with clean column names.
     Suitable for direct loading in R, Python, or Excel.

  2. pathway_results_llm.md
     Markdown document structured for input to a large language model (LLM).
     Contains a role instruction, biological context, column guide, top-50
     results table, seed gene list, reference set, and suggested questions.
     Paste into ChatGPT, Claude, or any LLM to discuss biological meaning.

Usage
-----
python pipeline/summarize_results.py \
    --scores        results/kf_chd/pathway_scores_standard.csv \
    --seeds         data/cohorts/chd/kf_chd_seeds.txt \
    --reference     data/cohorts/chd/kf_chd_cilia_reference.txt \
    --cohort-name   "KF-CHD" \
    --disease       "congenital heart disease" \
    --n-probands    697 \
    --outdir        results/kf_chd/

All arguments except --scores and --outdir are optional but strongly
recommended for a useful LLM output file.
"""

import argparse
import csv
import os
import sys
from datetime import date


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument('--scores',       required=True,
                   help='pathway_scores_standard.csv from score_pathways.py')
    p.add_argument('--seeds',        default=None,
                   help='Text file of seed gene CUIs (one per line)')
    p.add_argument('--edges-merged', default=None,
                   help='edges_merged.csv or edges_all_noncc.csv — used to compute '
                        'direct seed-pathway membership overlap (seed_members column). '
                        'Optional; if omitted, seed_members column is omitted.')
    p.add_argument('--nodes',        default=None,
                   help='nodes_clean_noncc.csv — used to map CUIs to gene symbols '
                        'in seed_members and contributing_seeds columns. Optional.')
    p.add_argument('--member-scores', default=None,
                   help='pathway_member_scores.tsv from score_pathways.py '
                        '--out-member-scores. Optional. Adds seed_member_scores column '
                        '(ppr_score_norm values parallel to seed_members, sorted by '
                        'PPR score descending).')
    p.add_argument('--influential-nodes', default=None,
                   help='pathway_influential_nodes.tsv from score_pathways.py '
                        '--out-influential-nodes. Optional. Adds influential_nodes_local '
                        '(z_within_pathway > 1.96) and influential_nodes_global '
                        '(z_global > 1.96) columns.')
    p.add_argument('--reference',    default=None,
                   help='Text file of reference pathway CUIs (one per line)')
    p.add_argument('--cohort-name',  default='BIFO run',
                   help='Short label, e.g. "KF-CHD"')
    p.add_argument('--disease',      default=None,
                   help='Disease or condition, e.g. "congenital heart disease"')
    p.add_argument('--n-probands',   type=int, default=None,
                   help='Number of probands in the cohort')
    p.add_argument('--n-resolved-seeds', type=int, default=None,
                   help='Number of seed genes resolved to graph CUIs (post-resolution). '
                        'If provided, shown in output instead of raw seed file count. '
                        'Use this to match the count used in the actual analysis.')
    p.add_argument('--outdir',       required=True,
                   help='Output directory')
    p.add_argument('--top-n-llm',    type=int, default=50,
                   help='Pathways in LLM table (default 50)')
    p.add_argument('--include-degenerate', action='store_true',
                   help='Include null-degenerate pathways in LLM table')
    return p.parse_args()


def load_scores(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def load_lines(path):
    if not path:
        return set()
    with open(path) as f:
        return {l.split()[0] for l in f if l.strip() and not l.startswith('#')}


def load_member_scores_file(path):
    """
    Load pathway_member_scores.tsv from score_pathways.py --out-member-scores.
    Returns {pathway_id: [(member_gene_name, ppr_score_norm), ...]} sorted by
    ppr_score_norm descending. Only rows for genes that appear in seed_members
    are used by summarize_results — the full file is kept for other uses.
    """
    if not path:
        return {}
    result = {}
    try:
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                pid  = row.get('pathway_id', '').strip()
                gene = row.get('member_gene_name', row.get('member_gene', '')).strip()
                norm = row.get('ppr_score_norm', '')
                if pid and gene and norm:
                    result.setdefault(pid, []).append(
                        (gene, float(norm) if norm else 0.0))
        # Sort each pathway's list by ppr_score_norm descending
        for pid in result:
            result[pid].sort(key=lambda x: -x[1])
    except FileNotFoundError:
        print(f"WARNING: member-scores file not found: {path}",
              file=__import__('sys').stderr)
    return result


def load_influential_nodes_file(path):
    """
    Load pathway_influential_nodes.tsv from score_pathways.py --out-influential-nodes.
    Returns:
      local_map  : {pathway_id: [node_name, ...]} where z_within_pathway > 1.96
      global_map : {pathway_id: [node_name, ...]} where z_global > 1.96
      z_stats    : {pathway_id: str} compact quantile summary of both z distributions,
                   format: "local:n=N;max=X;p75=X;p50=X;p25=X | global:n=N;max=X;p75=X;p50=X;p25=X"
    Both node lists sorted by their respective z-score descending.
    The z_stats string is always populated (even when local/global lists are empty),
    allowing users to distinguish "nothing passed threshold" from "file not loaded".
    """
    if not path:
        return {}, {}, {}
    local_map  = {}
    global_map = {}
    # Collect raw z-scores per pathway for quantile computation
    z_loc_raw  = {}   # {pathway_id: [z_within, ...]}
    z_glb_raw  = {}   # {pathway_id: [z_global, ...]}
    Z_THRESH = 1.96
    try:
        with open(path, newline='') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                pid   = row.get('pathway_id', '').strip()
                name  = row.get('node_name',  row.get('node_id', '')).strip()
                z_loc = row.get('z_within_pathway', '')
                z_glb = row.get('z_global', '')
                if not pid or not name:
                    continue
                try:
                    zl = float(z_loc)
                    z_loc_raw.setdefault(pid, []).append(zl)
                    if zl > Z_THRESH:
                        local_map.setdefault(pid, []).append((name, zl))
                except (ValueError, TypeError):
                    pass
                try:
                    zg = float(z_glb)
                    z_glb_raw.setdefault(pid, []).append(zg)
                    if zg > Z_THRESH:
                        global_map.setdefault(pid, []).append((name, zg))
                except (ValueError, TypeError):
                    pass
        for pid in local_map:
            local_map[pid].sort(key=lambda x: -x[1])
        for pid in global_map:
            global_map[pid].sort(key=lambda x: -x[1])
    except FileNotFoundError:
        print(f"WARNING: influential-nodes file not found: {path}",
              file=__import__('sys').stderr)
        return {}, {}, {}

    # Build quantile summary strings
    def _quantile(vals, q):
        """Simple quantile without numpy — linear interpolation."""
        s = sorted(vals)
        n = len(s)
        if n == 0:
            return float('nan')
        idx = q * (n - 1)
        lo, hi = int(idx), min(int(idx) + 1, n - 1)
        return s[lo] + (idx - lo) * (s[hi] - s[lo])

    def _fmt(v):
        return f"{v:.2f}"

    z_stats = {}
    all_pids = set(z_loc_raw) | set(z_glb_raw)
    for pid in all_pids:
        lv = z_loc_raw.get(pid, [])
        gv = z_glb_raw.get(pid, [])
        if lv:
            ls = (f"local:n={len(lv)};max={_fmt(max(lv))};"
                  f"p75={_fmt(_quantile(lv,0.75))};"
                  f"p50={_fmt(_quantile(lv,0.50))};"
                  f"p25={_fmt(_quantile(lv,0.25))}")
        else:
            ls = "local:n=0"
        if gv:
            gs = (f"global:n={len(gv)};max={_fmt(max(gv))};"
                  f"p75={_fmt(_quantile(gv,0.75))};"
                  f"p50={_fmt(_quantile(gv,0.50))};"
                  f"p25={_fmt(_quantile(gv,0.25))}")
        else:
            gs = "global:n=0"
        z_stats[pid] = f"{ls} | {gs}"

    return local_map, global_map, z_stats


def safe_float(val):
    try:
        f = float(val)
        return None if f != f else f  # f != f is True only for NaN
    except (TypeError, ValueError):
        return None


def safe_bool(val):
    if isinstance(val, str):
        return val.strip().lower() in ('true', '1', 'yes')
    return bool(val)


def build_cui_to_symbol(nodes_csv):
    """
    Build {CUI: symbol} from nodes CSV (plain or .gz). Returns empty dict if path is None.

    Also returns a secondary dict {any_node_id: CUI} mapping every node's 'label'
    field (which holds the source-vocabulary ID, e.g. 'HGNC:1097') back to its
    CUI (node_id / concept_id). This reverse map is used by build_membership_map
    to normalise gene member IDs from the edges file into CUIs regardless of which
    ID scheme the edges file uses for gene nodes (C-prefixed CUI vs HGNC: vs other).
    Call as: cui_to_symbol, nodeid_to_cui = build_cui_to_symbol(path)
    """
    if not nodes_csv:
        return {}, {}
    import csv as _csv, gzip
    lookup = {}        # CUI  -> symbol
    id_to_cui = {}     # any_node_id (label) -> CUI
    try:
        opener = gzip.open if nodes_csv.endswith('.gz') else open
        with opener(nodes_csv, 'rt', encoding='utf-8-sig') as f:
            for row in _csv.DictReader(f):
                cid = (row.get('node_id') or row.get('concept_id') or
                       row.get('id') or '').strip()
                sym = (row.get('name') or row.get('symbol') or
                       row.get('pref_term') or row.get('preferred_name') or '').strip()
                label = (row.get('label') or '').strip()
                # Skip header artifact rows
                if not cid or cid.startswith('-') or cid == 'CUI':
                    continue
                if cid and sym:
                    # Strip trailing ' gene' suffix from HGNC gene names
                    if sym.endswith(' gene'):
                        sym = sym[:-5].strip()
                    lookup[cid] = sym
                # Build reverse map: label (e.g. 'HGNC:1097') -> CUI
                if cid and label and label != cid:
                    id_to_cui[label] = cid
    except FileNotFoundError:
        print(f"WARNING: nodes file not found: {nodes_csv}", file=__import__('sys').stderr)
    return lookup, id_to_cui


MEMBERSHIP_PREDS_FWD = {
    'pathway_associated_with_gene',
    'has_signature_gene',
    'process_involves_gene',
}
MEMBERSHIP_PREDS_REV = {
    'inverse_pathway_associated_with_gene',
    'inverse_has_signature_gene',
    'gene_plays_role_in_process',
}


def build_membership_map(edges_csv, nodeid_to_cui=None, min_members=8, max_members=300):
    """
    Build {pathway_id: frozenset(gene_CUIs)} from edges_merged.csv.
    Only processes membership predicates; applies size filter.

    Gene member IDs from the edges file are normalised to CUIs using
    nodeid_to_cui (the reverse map from build_cui_to_symbol). This handles
    the case where MSIGDB/WikiPathways/Reactome membership edges store gene
    nodes under HGNC: or other vocabulary IDs rather than C-prefixed CUIs,
    which would otherwise prevent intersection with the seed CUI set.
    """
    import csv as _csv
    norm = nodeid_to_cui or {}
    pw_to_genes = {}
    try:
        import gzip as _gzip
        opener = _gzip.open if edges_csv.endswith('.gz') else open
        with opener(edges_csv, 'rt', encoding='utf-8-sig') as f:
            for row in _csv.DictReader(f):
                pred = row.get('predicate', '').strip()
                src  = row.get('source', row.get('subject', '')).strip()
                tgt  = row.get('target', row.get('object', '')).strip()
                # Strip trailing ' CUI' suffix from all node IDs
                src_clean = src[:-4].strip() if src.endswith(' CUI') else src
                tgt_clean = tgt[:-4].strip() if tgt.endswith(' CUI') else tgt
                # Normalise gene member IDs to CUIs via reverse lookup
                if pred in MEMBERSHIP_PREDS_FWD:
                    gene_id = norm.get(tgt_clean, tgt_clean)
                    pw_to_genes.setdefault(src_clean, set()).add(gene_id)
                elif pred in MEMBERSHIP_PREDS_REV:
                    gene_id = norm.get(src_clean, src_clean)
                    pw_to_genes.setdefault(tgt_clean, set()).add(gene_id)
    except FileNotFoundError:
        print(f"WARNING: edges file not found: {edges_csv}", file=__import__('sys').stderr)
        return {}
    return {pw: frozenset(genes)
            for pw, genes in pw_to_genes.items()
            if min_members <= len(genes) <= max_members}


def build_summary(rows, reference_ids, seed_ids=None, membership_map=None,
                  cui_to_symbol=None, member_scores_map=None,
                  influential_local_map=None, influential_global_map=None,
                  z_stats_map=None):
    seed_set = frozenset(seed_ids) if seed_ids else frozenset()
    sym = cui_to_symbol or {}
    ms_map   = member_scores_map   or {}
    inf_loc  = influential_local_map  or {}
    inf_glb  = influential_global_map or {}
    z_stats  = z_stats_map or {}

    def cuis_to_symbols(cui_str, delim='|'):
        """Convert pipe-separated CUI string to symbol string, falling back to CUI."""
        if not cui_str:
            return ''
        return ';'.join(sym.get(c.strip(), c.strip())
                        for c in cui_str.split(delim) if c.strip())

    def cui_set_to_symbols(cui_set):
        """Convert frozenset of CUIs to semicolon-separated symbol string."""
        return ';'.join(sorted(sym.get(c, c) for c in cui_set))

    for r in rows:
        r['_dn']  = safe_float(r.get('degree_norm', 0)) or 0.0
        r['_cal'] = safe_bool(r.get('null_calibrated', True))

    rows.sort(key=lambda r: r['_dn'], reverse=True)

    summary = []
    for rank, r in enumerate(rows, 1):
        nz  = safe_float(r.get('null_z'))
        eq  = safe_float(r.get('empirical_q'))
        mnz = safe_float(r.get('member_mean_null_z'))
        mq  = safe_float(r.get('member_mean_q'))
        cid = r.get('concept_id', '')

        # seed_members: CUIs resolved to symbols, sorted by PPR score if available
        raw_seed_member_cuis = (
            membership_map.get(cid, frozenset()) & seed_set
            if (membership_map is not None and seed_set) else frozenset()
        )
        # Build ordered seed member list: by ppr_score_norm desc if ms_map available,
        # otherwise alphabetically
        if ms_map and cid in ms_map:
            # ms_map[cid] is [(gene_name, ppr_score_norm), ...] sorted desc
            # Convert raw_seed_member_cuis to symbols for matching
            seed_sym_set = {sym.get(c, c) for c in raw_seed_member_cuis}
            ordered_seed_members = [
                gene for gene, _ in ms_map[cid]
                if gene in seed_sym_set
            ]
            # Append any seed members not in ms_map (shouldn't happen, but safe)
            covered = set(ordered_seed_members)
            for c in sorted(raw_seed_member_cuis):
                s = sym.get(c, c)
                if s not in covered:
                    ordered_seed_members.append(s)
            seed_members_str = ';'.join(ordered_seed_members)
            # seed_member_scores: parallel ppr_score_norm values
            ms_lookup = {gene: norm for gene, norm in ms_map[cid]}
            seed_member_scores_str = ';'.join(
                f"{ms_lookup.get(g, ''):.4f}" if ms_lookup.get(g) is not None else ''
                for g in ordered_seed_members
            )
        else:
            seed_members_str = cui_set_to_symbols(raw_seed_member_cuis)
            seed_member_scores_str = ''

        # influential nodes: z-thresholded, semicolon-separated node names
        inf_local_str = ';'.join(
            name for name, _ in inf_loc.get(cid, [])
        )
        inf_global_str = ';'.join(
            name for name, _ in inf_glb.get(cid, [])
        )

        summary.append({
            'rank':                    rank,
            'pathway_name':            r.get('name', cid),
            'pathway_id':              cid,
            'source':                  r.get('sab', ''),
            'n_members':               r.get('member_gene_count') or r.get('n_members') or r.get('degree') or '0',
            'degree_norm':             f"{r['_dn']:.6e}",
            'null_calibrated':         str(r['_cal']),
            'null_z':                  f"{nz:.3f}"  if nz  is not None else 'NaN',
            'empirical_q':             f"{eq:.4f}"  if eq  is not None else 'NaN',
            'member_mean_null_z':      f"{mnz:.3f}" if mnz is not None else 'NaN',
            'member_mean_q':           f"{mq:.4f}"  if mq  is not None else 'NaN',
            'in_reference':            str(cid in reference_ids) if reference_ids else 'NA',
            'contributing_seeds':      cuis_to_symbols(r.get('contributing_seeds', '')),
            'seed_members':            seed_members_str,
            'seed_member_scores':      seed_member_scores_str,
            'influential_nodes_local': inf_local_str,
            'influential_nodes_global': inf_global_str,
            'neighbor_z_summary':      z_stats.get(cid, ''),
            '_cal':                    r['_cal'],
        })
    return summary


TSV_COLS = [
    'rank', 'pathway_name', 'pathway_id', 'source', 'n_members',
    'degree_norm', 'null_calibrated', 'null_z', 'empirical_q',
    'member_mean_null_z', 'member_mean_q', 'in_reference',
    'contributing_seeds', 'seed_members', 'seed_member_scores',
    'influential_nodes_local', 'influential_nodes_global',
    'neighbor_z_summary',
]


def write_tsv(summary, outpath):
    with open(outpath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=TSV_COLS, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        w.writerows(summary)
    print(f"Saved: {outpath}  ({len(summary)} pathways)")


LLM_COLS = [
    'rank', 'pathway_name', 'source', 'n_members',
    'degree_norm', 'null_z', 'empirical_q', 'in_reference',
    'contributing_seeds', 'seed_members', 'seed_member_scores',
    'influential_nodes_local', 'influential_nodes_global',
    'neighbor_z_summary',
]


def md_table(rows, cols):
    widths = {c: len(c) for c in cols}
    for r in rows:
        for c in cols:
            widths[c] = max(widths[c], len(str(r.get(c, ''))))
    def row_str(r):
        return '| ' + ' | '.join(str(r.get(c,'')).ljust(widths[c]) for c in cols) + ' |'
    hdr = '| ' + ' | '.join(c.ljust(widths[c]) for c in cols) + ' |'
    sep = '| ' + ' | '.join('-'*widths[c] for c in cols) + ' |'
    return '\n'.join([hdr, sep] + [row_str(r) for r in rows])


def write_llm(summary, seed_ids, reference_ids, cohort_name, disease,
              n_probands, scores_path, top_n, include_degenerate, outpath,
              n_resolved_seeds=None):

    n_seeds      = len(seed_ids)
    n_display    = n_resolved_seeds if n_resolved_seeds is not None else n_seeds
    seeds_note   = (f' ({n_seeds} in input file; {n_display} resolved to graph CUIs)'
                    if n_resolved_seeds is not None else '')
    n_total      = len(summary)
    n_cal        = sum(1 for r in summary if r['_cal'])
    n_sig        = sum(1 for r in summary
                       if safe_float(r['empirical_q']) is not None
                       and safe_float(r['empirical_q']) < 0.05)
    disease_str  = disease or 'the condition of interest'
    proband_str  = f"{n_probands} probands" if n_probands else "a patient cohort"
    today        = date.today().isoformat()

    top_rows = (summary if include_degenerate
                else [r for r in summary if r['_cal']])[:top_n]

    top_nz = sorted(
        [r for r in summary if safe_float(r['null_z']) is not None],
        key=lambda r: safe_float(r['null_z']), reverse=True
    )[:5]

    sig_rows = [r for r in summary
                if safe_float(r['empirical_q']) is not None
                and safe_float(r['empirical_q']) < 0.05]

    seed_list = sorted(seed_ids)
    if len(seed_list) > 200:
        seed_display = ', '.join(seed_list[:200]) + \
                       f'\n... and {len(seed_list)-200} more (see --seeds file).'
    else:
        seed_display = ', '.join(seed_list) if seed_list else 'Not provided.'

    ref_display = (', '.join(sorted(reference_ids))
                   if reference_ids else 'No reference set provided.')

    L = []

    L.append('## Instructions for AI\n')
    L.append(
        'You are helping a scientist interpret the **biological meaning** of pathway '
        'analysis results from a computational tool called BIFO (Biological Information '
        'Flow Ontology). Your role is to:\n\n'
        '- Help the user understand which biological processes are enriched in their data\n'
        '- Explain what the pathway names mean in plain language\n'
        '- Connect the results to the disease, condition, or experimental system the user is studying\n'
        '- Suggest follow-up questions, experiments, or analyses that might be warranted\n'
        '- Answer questions about the input data (what genes were used, what the analysis did)\n\n'
        'This session is intended for biological interpretation of the results. '
        'Unless the user asks you to evaluate the method, please focus this conversation '
        'on the biology: what the pathways mean, how they relate to the disease or '
        'experimental system, and what follow-up questions or experiments might be '
        'worth pursuing. If the user wants to discuss the method, answer their questions, '
        'then return focus to the biology.\n\n'
        '*Note: these outputs are intended for exploratory biological interpretation. '
        'They are not part of any quantitative evaluation and are not used for '
        'scoring, ranking, or statistical inference.*\n\n'
        'The user may be a biologist, clinician, or computational scientist. '
        'Do not assume prior knowledge of graph algorithms or bioinformatics. '
        'The user may have questions specific to their own experimental system, '
        'organism, tissue type, or clinical context, not necessarily the disease listed below.\n'
    )
    L.append('---\n')

    L.append(f'## Analysis Summary: {cohort_name}\n\n'
             f'- **Cohort:** {cohort_name}\n'
             f'- **Disease / condition:** {disease_str}\n'
             f'- **Cohort size:** {proband_str}\n'
             f'- **Input genes (seeds):** {n_display}{seeds_note}\n'
             f'- **Pathways scored:** {n_total} total; {n_cal} with valid null distribution\n'
             f'- **Significant pathways (q < 0.05):** {n_sig}\n'
             f'- **Analysis date:** {today}\n'
             f'- **Source file:** {os.path.basename(scores_path)}\n')
    L.append('---\n')

    L.append(
        '## What This Analysis Does\n\n'
        'BIFO propagates signal from a set of input genes through a large biomedical '
        'knowledge graph connecting genes, proteins, pathways, diseases, and biological '
        'processes. Only biologically meaningful relationships are used; statistical '
        'associations and text-mining edges are excluded.\n\n'
        'Signal flows from the input genes outward and accumulates at pathway nodes. '
        'Pathways strongly connected to the input genes through biological relationships '
        'receive higher scores. Scores are then compared to a statistical null model '
        '(1,000 random rewirings of gene-pathway membership) to identify pathways that '
        'score higher than expected by chance given the graph topology.\n\n'
        'This approach is particularly useful when the input gene list is large or '
        'heterogeneous because it concentrates distributed biological signal rather than '
        'requiring strong direct overlap between input genes and pathway members.\n'
    )
    L.append('---\n')

    L.append(
        '## How to Read the Results Table\n\n'
        '| Column | What it means |\n'
        '|--------|---------------|\n'
        '| **rank** | Position by BIFO score (1 = highest). Lower rank = more signal from input genes. |\n'
        '| **pathway_name** | Name of the biological pathway or process. |\n'
        '| **source** | Database the pathway comes from (MSIGDB, WikiPathways, Reactome, GO). |\n'
        '| **n_members** | Number of genes in the pathway. |\n'
        '| **degree_norm** | BIFO score: propagated signal at the pathway node, normalised by pathway size. |\n'
        '| **null_z** | Standard deviations above expected score vs. 1,000 random rewirings. Higher = more specific enrichment. NaN = test not valid. |\n'
        '| **empirical_q** | BH-corrected p-value. q < 0.05 = statistically significant. NaN = not tested. |\n'
        '| **in_reference** | TRUE if in a pre-specified biologically relevant reference set. |\n'
        '| **contributing_seeds** | Seed genes whose PPR signal flowed through the graph and contributed to this pathway\'s score. Includes genes connected through mechanistic edges, not just direct pathway members. |\n'
        '| **seed_members** | Seed genes that are direct annotated members of this pathway in the source database. Sorted by PPR score descending when scoring data is available. This is what Fisher enrichment uses. |\n'
        '| **seed_member_scores** | PPR score (normalised 0-1 within pathway, 1.0 = top gene) for each gene in seed_members, in the same order. Empty if --member-scores not provided. |\n'
        '| **influential_nodes_local** | Nodes in the conditioned graph connected to this pathway whose PPR score exceeds z > 1.96 within this pathway\'s neighbor distribution. Includes non-seed, non-member nodes that drove the pathway\'s score through graph propagation. |\n'
        '| **influential_nodes_global** | Same nodes filtered by z > 1.96 vs all HGNC gene nodes in the full graph. A stricter threshold: these nodes are globally high-scoring, not just locally prominent for this pathway. |\n'
        '| **neighbor_z_summary** | Quantile summary of z-score distributions across all neighbor nodes for this pathway, for both local (within-pathway) and global (all genes) z-scores. Format: local:n=N;max=X;p75=X;p50=X;p25=X | global:n=N;max=X;p75=X;p50=X;p25=X. Use this to interpret empty influential_nodes columns: if max z is below 1.96, the threshold was not met rather than the data being missing. |\n\n'
        '**rank and null_z are complementary, not interchangeable.** '
        'Rank reflects absolute propagated signal; null_z reflects whether signal is '
        'specifically attributable to this pathway\'s gene membership rather than generic '
        'graph topology. null_z reflects separation from a graph-specific null and is '
        'not a directly comparable effect size across different analyses.\n'
    )
    L.append('---\n')

    cal_note = '' if include_degenerate else ' (well-calibrated pathways only)'
    L.append(f'## Top {top_n} Pathways by BIFO Score{cal_note}\n\n'
             f'*Full results in pathway_results_summary.tsv*\n\n'
             + md_table(top_rows, LLM_COLS))
    L.append('\n---\n')

    L.append(f'## Key Findings\n\n**Significant pathways (q < 0.05):** {len(sig_rows)}\n')
    sig_shown = set()
    if sig_rows:
        for r in sig_rows[:10]:
            L.append(f"- {r['pathway_name']} "
                     f"(rank {r['rank']}, null_z = {r['null_z']}, q = {r['empirical_q']})")
            sig_shown.add(r['pathway_id'])
        L.append('')
    L.append('\n**Top 5 by null_z (strongest statistical enrichment):**\n')
    for r in top_nz:
        if r['pathway_id'] not in sig_shown:
            L.append(f"- {r['pathway_name']} "
                     f"(rank {r['rank']}, null_z = {r['null_z']}, q = {r['empirical_q']})")
    L.append('\n---\n')

    L.append(
        '## Suggested Questions to Ask\n\n'
        f'- "What biological processes are enriched in {disease_str}?"\n'
        '- "What do the top-ranked pathways have in common biologically?"\n'
        '- "Which pathways are both highly ranked and statistically significant?"\n'
        '- "Are any of these pathways relevant to [your tissue / cell type / organism]?"\n'
        '- "What genes from my input list are members of the top pathways?"\n'
        '- "How might these pathway findings inform follow-up experiments?"\n'
        '- "What is the biological function of [pathway name]?"\n'
        '- "What were the input genes used to generate these results?"\n'
        '- "What database did pathway [X] come from, and how was it defined?"\n'
        '- "Are any of these pathways relevant to [specific disease or condition]?"\n'
    )
    L.append('---\n')

    L.append(
        f'## Input Data Context\n\n'
        f'### Seed genes (input to BIFO) — {n_display} resolved{seeds_note}\n\n'
        'These genes were used as the starting point for signal propagation. '
        'They are typically genes carrying rare variants in the cohort, or genes '
        'of interest in the experimental system.\n\n'
        + seed_display + '\n\n'
        '### Reference pathway set\n\n'
        'These pathways were pre-specified as biologically relevant before the '
        'analysis was run. They are marked TRUE in the in_reference column.\n\n'
        + ref_display + '\n\n'
        '### Analysis parameters\n\n'
        f'- Scores file: {os.path.basename(scores_path)}\n'
        f'- Total pathways scored: {n_total}\n'
        f'- Pathways with valid null distribution: {n_cal}\n'
        '- Null model: degree-preserving membership rewiring, N = 1,000 permutations\n'
        '- Multiple testing correction: Benjamini-Hochberg\n'
        '- Pathway size filter: minimum 8 members, maximum 300 members\n'
        '- Pathway sources: MSigDB canonical collections, WikiPathways, Reactome, GO\n\n'
        f'---\n*Generated by summarize_results.py on {today}*\n'
    )

    with open(outpath, 'w') as f:
        f.write('\n'.join(L))
    print(f"Saved: {outpath}  (top {top_n} pathways in LLM table)")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading: {args.scores}")
    rows = load_scores(args.scores)
    if not rows:
        print("ERROR: no rows in scores file", file=sys.stderr)
        sys.exit(1)

    seed_ids      = load_lines(args.seeds)
    reference_ids = load_lines(args.reference)
    print(f"  {len(rows)} pathways, {len(seed_ids)} seeds, "
          f"{len(reference_ids)} reference pathways")

    cui_to_symbol, nodeid_to_cui = build_cui_to_symbol(args.nodes)
    if cui_to_symbol:
        print(f"  {len(cui_to_symbol):,} CUI->symbol mappings loaded")
    if nodeid_to_cui:
        print(f"  {len(nodeid_to_cui):,} node-id->CUI reverse mappings loaded")

    membership_map = None
    if args.edges_merged:
        print(f"Building membership map from: {args.edges_merged}")
        membership_map = build_membership_map(args.edges_merged,
                                              nodeid_to_cui=nodeid_to_cui)
        print(f"  {len(membership_map):,} pathways in membership map")

    member_scores_map = None
    if args.member_scores:
        print(f"Loading member scores from: {args.member_scores}")
        member_scores_map = load_member_scores_file(args.member_scores)
        print(f"  {len(member_scores_map):,} pathways with member score data")

    influential_local_map  = None
    influential_global_map = None
    z_stats_map            = None
    if args.influential_nodes:
        print(f"Loading influential nodes from: {args.influential_nodes}")
        influential_local_map, influential_global_map, z_stats_map = \
            load_influential_nodes_file(args.influential_nodes)
        print(f"  {len(influential_local_map):,} pathways with local z>1.96 nodes, "
              f"{len(influential_global_map):,} with global z>1.96 nodes")

    summary = build_summary(rows, reference_ids,
                            seed_ids=seed_ids,
                            membership_map=membership_map,
                            cui_to_symbol=cui_to_symbol,
                            member_scores_map=member_scores_map,
                            influential_local_map=influential_local_map,
                            influential_global_map=influential_global_map,
                            z_stats_map=z_stats_map)

    write_tsv(summary, os.path.join(args.outdir, 'pathway_results_summary.tsv'))
    write_llm(
        summary            = summary,
        seed_ids           = seed_ids,
        reference_ids      = reference_ids,
        cohort_name        = args.cohort_name,
        disease            = args.disease,
        n_probands         = args.n_probands,
        scores_path        = args.scores,
        top_n              = args.top_n_llm,
        include_degenerate = args.include_degenerate,
        n_resolved_seeds   = args.n_resolved_seeds,
        outpath            = os.path.join(args.outdir, 'pathway_results_llm.md'),
    )
    print("Done.")


if __name__ == '__main__':
    main()
