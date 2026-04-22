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
    p.add_argument('--reference',    default=None,
                   help='Text file of reference pathway CUIs (one per line)')
    p.add_argument('--cohort-name',  default='BIFO run',
                   help='Short label, e.g. "KF-CHD"')
    p.add_argument('--disease',      default=None,
                   help='Disease or condition, e.g. "congenital heart disease"')
    p.add_argument('--n-probands',   type=int, default=None,
                   help='Number of probands in the cohort')
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
        return {l.strip() for l in f if l.strip() and not l.startswith('#')}


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


def build_summary(rows, reference_ids):
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
        summary.append({
            'rank':               rank,
            'pathway_name':       r.get('name', cid),
            'pathway_id':         cid,
            'source':             r.get('sab', ''),
            'n_members':          r.get('member_gene_count') or r.get('n_members') or r.get('degree') or '0',
            'degree_norm':        f"{r['_dn']:.6e}",
            'null_calibrated':    str(r['_cal']),
            'null_z':             f"{nz:.3f}"  if nz  is not None else 'NaN',
            'empirical_q':        f"{eq:.4f}"  if eq  is not None else 'NaN',
            'member_mean_null_z': f"{mnz:.3f}" if mnz is not None else 'NaN',
            'member_mean_q':      f"{mq:.4f}"  if mq  is not None else 'NaN',
            'in_reference':       str(cid in reference_ids) if reference_ids else 'NA',
            '_cal':               r['_cal'],
        })
    return summary


TSV_COLS = [
    'rank', 'pathway_name', 'pathway_id', 'source', 'n_members',
    'degree_norm', 'null_calibrated', 'null_z', 'empirical_q',
    'member_mean_null_z', 'member_mean_q', 'in_reference',
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
              n_probands, scores_path, top_n, include_degenerate, outpath):

    n_seeds      = len(seed_ids)
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
             f'- **Input genes (seeds):** {n_seeds}\n'
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
        '| **in_reference** | TRUE if in a pre-specified biologically relevant reference set. |\n\n'
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
        f'### Seed genes (input to BIFO) — {n_seeds} total\n\n'
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

    summary = build_summary(rows, reference_ids)

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
        outpath            = os.path.join(args.outdir, 'pathway_results_llm.md'),
    )
    print("Done.")


if __name__ == '__main__':
    main()
