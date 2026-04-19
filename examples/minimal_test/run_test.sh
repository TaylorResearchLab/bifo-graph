#!/usr/bin/env bash
# =============================================================================
# run_test.sh — Minimal BIFO pipeline test run
# =============================================================================
# Prerequisites:
#   - nodes.csv and edges_raw.csv exported from Neo4j (see neo4j_export.cypher)
#   - seed_nodes.txt and heldout_nodes.txt populated with real CUIs
#   - bifo_ddkg_mapping.yaml v0.7.1 present
#   - bifo_conditioning.py and score_pathways.py on PATH or in this directory
#
# Usage:
#   bash run_test.sh
#
# All outputs written to ./test_output/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/test_output"
mkdir -p "${OUT_DIR}"

# --- Paths (edit if your files are elsewhere) ---
NODES="${SCRIPT_DIR}/nodes.csv"
EDGES="${SCRIPT_DIR}/edges_raw.csv"
MEMBERSHIP_EDGES="${SCRIPT_DIR}/pathway_membership_edges.csv"
EDGES_MERGED="${OUT_DIR}/edges_merged.csv"
MAPPING="${SCRIPT_DIR}/bifo_ddkg_mapping.yaml"
SEEDS="${SCRIPT_DIR}/seed_nodes.txt"
HELDOUT="${SCRIPT_DIR}/heldout_nodes.txt"
CHD_REF="${SCRIPT_DIR}/chd_pathway_reference.txt"
BIFO_COND="${SCRIPT_DIR}/bifo_conditioning.py"
SCORE_PW="${SCRIPT_DIR}/score_pathways.py"

# --- Output stem ---
OUT_JSON="${OUT_DIR}/results.json"

echo "============================================================"
echo "BIFO Minimal Test Run"
echo "============================================================"
echo "Nodes:      ${NODES}"
echo "Edges:      ${EDGES}"
echo "Mapping:    ${MAPPING}"
echo "Seeds:      $(grep -v '^#' ${SEEDS} | grep -c '[^[:space:]]' || echo '?') genes"
echo "Held-out:   $(grep -v '^#' ${HELDOUT} | grep -c '[^[:space:]]' || echo '?') genes"
echo "Output dir: ${OUT_DIR}"
echo ""

# --- Preflight: check nodes.csv before scoring ---
echo "[0/3] Preflight check on nodes.csv..."
python3 "${SCORE_PW}" \
    --nodes "${NODES}" \
    --preflight \
    | python3 -c "
import json, sys
report = json.load(sys.stdin)
print(f\"  Concepts:      {report.get('n_concepts', '?')}\")
print(f\"  Multi-SAB:     {report.get('n_multi_sab_concepts', '?')}\")
print(f\"  Pathway SABs:  {report.get('pathway_sab_counts', {})}\")
print(f\"  Advice:        {report.get('recommendation', '?')}\")
"
echo ""

# --- Merge edges: combine seed<->hop1 edges with targeted membership edges ---
echo "[1/3] Merging edges_raw.csv + pathway_membership_edges.csv..."
if [ -f "${MEMBERSHIP_EDGES}" ]; then
    cp "${EDGES}" "${EDGES_MERGED}"
    # Extract only source,target,predicate columns (columns 1,2,3) from membership
    # file before appending — Neo4j exports may include extra metadata columns
    python3 - <<PYEOF
import csv, sys
with open("${MEMBERSHIP_EDGES}", encoding="utf-8-sig") as fin, open("${EDGES_MERGED}", 'a') as fout:
    reader = csv.DictReader(fin)
    # Normalise column names to lowercase for robustness
    fieldnames = [f.lower().strip() for f in (reader.fieldnames or [])]
    # Map to expected names
    src_col  = next((f for f in fieldnames if f in ('source','src','a','start')), fieldnames[0] if fieldnames else None)
    tgt_col  = next((f for f in fieldnames if f in ('target','tgt','b','end')),   fieldnames[1] if len(fieldnames)>1 else None)
    pred_col = next((f for f in fieldnames if f in ('predicate','type','relation','edge_type')), fieldnames[2] if len(fieldnames)>2 else None)
    if not (src_col and tgt_col and pred_col):
        print(f"ERROR: could not identify source/target/predicate columns in {'"'"'${MEMBERSHIP_EDGES}'"'"'}", file=sys.stderr)
        print(f"  Found columns: {fieldnames}", file=sys.stderr)
        sys.exit(1)
    # Ensure edges_raw.csv ends with a newline before appending
    fout.write("\n")
    for row in reader:
        cols = [f.lower().strip() for f in row.keys()]
        row_lower = {f.lower().strip(): v for f, v in row.items()}
        fout.write(f"{row_lower[src_col]},{row_lower[tgt_col]},{row_lower[pred_col]}\n")
PYEOF
    EDGE_COUNT=$(wc -l < "${EDGES_MERGED}")
    echo "  Merged edge file: ${EDGE_COUNT} lines -> ${EDGES_MERGED}"
else
    echo "  WARNING: pathway_membership_edges.csv not found -- using edges_raw.csv only"
    echo "  Run Step 4 in neo4j_export.cypher to generate membership edges"
    cp "${EDGES}" "${EDGES_MERGED}"
fi
echo ""

# --- Step 2: Conditioning + propagation (full merged graph) ---
echo "[2/4] Running bifo_conditioning.py — FULL graph (edges_raw + membership edges)..."
python3 "${BIFO_COND}" \
    --nodes   "${NODES}" \
    --edges   "${EDGES_MERGED}" \
    --mapping "${MAPPING}" \
    --seed-nodes    "${SEEDS}" \
    --heldout-nodes "${HELDOUT}" \
    --out-json "${OUT_DIR}/results.json" \
    --alpha 0.5

echo ""
echo "  Full graph conditioning outputs:"
ls -lh "${OUT_DIR}"/results_* 2>/dev/null | grep -v ablation

python3 - <<PYEOF
import json
with open("${OUT_DIR}/results.json") as f:
    r = json.load(f)
print(f"  Graph stats: {r.get('graph_stats', {})}")
for arm in ['raw', 'conditioned']:
    m = r.get(arm, {})
    print(f"  {arm:25s}: auroc={m.get('auroc','nan'):.3f}  "
          f"auprc={m.get('auprc','nan'):.3f}  "
          f"entropy={m.get('entropy','nan'):.3f}")
PYEOF
echo ""

# --- Step 2b: True propagation ablation — conditioning on edges_raw ONLY ---
echo "[2b/4] Running bifo_conditioning.py — ABLATION graph (edges_raw only)..."
echo "  NOTE: This is the true propagation ablation."
echo "  PPR is computed without membership edges."
echo "  Pathway scoring (step 3b) will use the full membership map."
python3 "${BIFO_COND}" \
    --nodes   "${NODES}" \
    --edges   "${EDGES}" \
    --mapping "${MAPPING}" \
    --seed-nodes    "${SEEDS}" \
    --heldout-nodes "${HELDOUT}" \
    --out-json "${OUT_DIR}/results_ablation.json" \
    --alpha 0.5

echo ""
echo "  Ablation conditioning outputs:"
ls -lh "${OUT_DIR}"/results_ablation_* 2>/dev/null

python3 - <<PYEOF
import json
with open("${OUT_DIR}/results_ablation.json") as f:
    r = json.load(f)
for arm in ['raw', 'conditioned']:
    m = r.get(arm, {})
    print(f"  {arm:25s}: auroc={m.get('auroc','nan'):.3f}  "
          f"auprc={m.get('auprc','nan'):.3f}  "
          f"entropy={m.get('entropy','nan'):.3f}")
PYEOF
echo ""

# --- Step 2c: Mechanistic-only arm ---
echo "[2c/4] Running bifo_conditioning.py — MECHANISTIC-ONLY..."
echo "  PPR uses only classification=='mechanistic' edges:"
echo "  Signal Transduction, Transcription, Translation, Protein Interaction."
echo "  Pathway Contribution (~85%), Observational, Perturbational all excluded."
python3 "${BIFO_COND}" \
    --nodes   "${NODES}" \
    --edges   "${EDGES_MERGED}" \
    --mapping "${MAPPING}" \
    --seed-nodes    "${SEEDS}" \
    --heldout-nodes "${HELDOUT}" \
    --out-json "${OUT_DIR}/results_mech.json" \
    --alpha 0.5 \
    --mechanistic-only

echo ""
echo "  Mechanistic-only outputs:"
ls -lh "${OUT_DIR}"/results_mech_* 2>/dev/null

python3 - <<PYEOF
import json
with open("${OUT_DIR}/results_mech.json") as f:
    r = json.load(f)
gs = r.get('graph_stats', {})
print(f"  Propagating edges: {gs.get('conditioned_propagating_edges', '?')}")
fd = gs.get('flow_class_distribution', {})
for k, v in sorted(fd.items(), key=lambda x: -x[1]):
    print(f"    {k:35s} {v:6,}")
for arm in ['raw', 'conditioned']:
    m = r.get(arm, {})
    print(f"  {arm:25s}: auroc={m.get('auroc','nan'):.3f}  "
          f"auprc={m.get('auprc','nan'):.3f}  "
          f"entropy={m.get('entropy','nan'):.3f}")
PYEOF
echo ""

# --- Step 3: Pathway scoring — three arms ---
echo "[3/4] Running score_pathways.py (Analysis B)..."
echo ""

# Arm A: full graph PPR + full membership (current best)
echo "  [3a] Full graph PPR + full membership (min=8, max=300)..."
python3 "${SCORE_PW}" \
    --nodes              "${NODES}" \
    --edges-raw          "${EDGES_MERGED}" \
    --edges-conditioned  "${OUT_DIR}/results_kept_edges.csv" \
    --scores-cond        "${OUT_DIR}/results_scores_cond.npy" \
    --scores-raw         "${OUT_DIR}/results_scores_raw.npy" \
    --node-index         "${OUT_DIR}/results_node_index.json" \
    --seed-nodes         "${SEEDS}" \
    --chd-pathways       "${CHD_REF}" \
    --min-members        8 \
    --max-members        300 \
    --out-csv            "${OUT_DIR}/pathway_scores_full.csv" \
    --out-json           "${OUT_DIR}/pathway_metrics_full.json" \
    --score-variant      degree_norm
echo ""

# Baseline comparison for CHD curated benchmark (Analysis 4)
echo "  [3a-baseline] CHD baseline enrichment comparison..."
python3 "${SCRIPT_DIR}/baseline_enrichment.py" \
    --edges-merged  "${EDGES_MERGED}" \
    --node-index    "${OUT_DIR}/results_node_index.json" \
    --scores-raw    "${OUT_DIR}/results_scores_raw.npy" \
    --scores-cond   "${OUT_DIR}/results_scores_cond.npy" \
    --bifo-scores   "${OUT_DIR}/pathway_scores_full.csv" \
    --chd-pathways  "${CHD_REF}" \
    --seed-nodes    "${SEEDS}" \
    --out-csv       "${OUT_DIR}/baseline_comparison.csv" \
    --out-json      "${OUT_DIR}/baseline_comparison.json" \
    --min-members   8 --max-members 300
echo ""
# This tests whether propagation signal without membership edges is sufficient
echo "  [3b] TRUE ABLATION: edges_raw PPR + full membership (min=8, max=300)..."
python3 "${SCORE_PW}" \
    --nodes              "${NODES}" \
    --edges-raw          "${EDGES_MERGED}" \
    --edges-conditioned  "${OUT_DIR}/results_ablation_kept_edges.csv" \
    --scores-cond        "${OUT_DIR}/results_ablation_scores_cond.npy" \
    --scores-raw         "${OUT_DIR}/results_ablation_scores_raw.npy" \
    --node-index         "${OUT_DIR}/results_ablation_node_index.json" \
    --seed-nodes         "${SEEDS}" \
    --chd-pathways       "${CHD_REF}" \
    --min-members        8 \
    --max-members        300 \
    --out-csv            "${OUT_DIR}/pathway_scores_ablation.csv" \
    --out-json           "${OUT_DIR}/pathway_metrics_ablation.json" \
    --score-variant      degree_norm
echo ""

# Arm C: mechanistic-only PPR + full membership
echo "  [3c] MECHANISTIC-ONLY PPR + full membership (min=8, max=300)..."
python3 "${SCORE_PW}" \
    --nodes              "${NODES}" \
    --edges-raw          "${EDGES_MERGED}" \
    --edges-conditioned  "${OUT_DIR}/results_mech_kept_edges.csv" \
    --scores-cond        "${OUT_DIR}/results_mech_scores_cond.npy" \
    --scores-raw         "${OUT_DIR}/results_mech_scores_raw.npy" \
    --node-index         "${OUT_DIR}/results_mech_node_index.json" \
    --seed-nodes         "${SEEDS}" \
    --chd-pathways       "${CHD_REF}" \
    --min-members        8 \
    --max-members        300 \
    --out-csv            "${OUT_DIR}/pathway_scores_mech.csv" \
    --out-json           "${OUT_DIR}/pathway_metrics_mech.json" \
    --score-variant      degree_norm
echo ""

# Summary — all three arms
echo "  [summary] Three-arm comparison..."
python3 - <<PYEOF
import json

def show(label, path):
    try:
        m = json.load(open(path))['metrics']
        n   = int(m.get('n_pathways_scored', 0))
        chd = int(m.get('n_chd_in_reference', 0))
        p10 = m.get('top10_precision', 0)
        e10 = m.get('top10_enrichment_ratio', float('nan'))
        ri  = m.get('rank_improvement_cond_vs_raw', float('nan'))
        mr  = m.get('mean_rank_chd_cond', float('nan'))
        print(f"  {label:42s}  n={n:4d}  CHD={chd:2d}  top10={p10:.2f}  "
              f"enrich={e10:.1f}x  mean_rank={mr:.0f}  rank_imp={ri:+.1f}")
    except Exception as e:
        print(f"  {label:42s}  ERROR: {e}")

show("Full graph PPR + full membership",      "${OUT_DIR}/pathway_metrics_full.json")
show("Ablation (edges_raw) + full membership","${OUT_DIR}/pathway_metrics_ablation.json")
show("Mechanistic-only PPR + full membership","${OUT_DIR}/pathway_metrics_mech.json")
PYEOF

# Baseline summary
echo ""
echo "  [baseline summary] CHD baseline comparison..."
python3 - <<PYEOF
import json

def show_baseline(label, path):
    try:
        data = json.load(open(path))
        for m in data.get('methods', []):
            method = m.get('method','?')
            p10  = m.get('precision_at_10', 0)
            ndcg = m.get('ndcg_at_10', 0)
            ap   = m.get('average_precision', 0)
            mr   = m.get('mean_ref_rank', float('nan'))
            print(f"  {label:15s}  {method:25s}  P@10={p10:.3f}  NDCG@10={ndcg:.4f}  AP={ap:.4f}  MeanRank={mr:.1f}")
    except Exception as e:
        print(f"  {label:15s}  ERROR: {e}")

show_baseline("CHD curated",  "${OUT_DIR}/baseline_comparison.json")
PYEOF

echo ""
echo "============================================================"
echo "Run complete."
echo "  Gene recovery (full):     ${OUT_DIR}/results.json"
echo "  Gene recovery (ablation): ${OUT_DIR}/results_ablation.json"
echo "  Pathway scores (full):    ${OUT_DIR}/pathway_scores_full.csv"
echo "  Pathway scores (ablation):${OUT_DIR}/pathway_scores_ablation.csv"
echo "  Ablation isolates propagation contribution:"
echo "  both arms use identical membership maps for scoring."
echo "============================================================"

# =============================================================================
# C4: Reactome pathway-split controls
# =============================================================================
# Uses pathway gene members as seeds (70%) and heldout (30%).
# Tests whether BIFO recovers the source pathway and related pathways
# without relying on hand-curated disease gene lists.
# Two splits: Notch signaling (cardiac-relevant) and MAPK (cancer-relevant).
# =============================================================================

run_c4() {
    local LABEL="$1"
    local SEEDS_C4="$2"
    local HELDOUT_C4="$3"
    local REF_C4="$4"
    local OUT_C4="${OUT_DIR}/c4_${LABEL}"

    mkdir -p "${OUT_C4}"
    echo ""
    echo "  [C4/${LABEL}] Conditioning..."
    python3 "${BIFO_COND}" \
        --nodes   "${NODES}" \
        --edges   "${EDGES_MERGED}" \
        --mapping "${MAPPING}" \
        --seed-nodes    "${SEEDS_C4}" \
        --heldout-nodes "${HELDOUT_C4}" \
        --out-json "${OUT_C4}/results.json" \
        --alpha 0.5

    echo "  [C4/${LABEL}] Pathway scoring..."
    python3 "${SCORE_PW}" \
        --nodes              "${NODES}" \
        --edges-raw          "${EDGES_MERGED}" \
        --edges-conditioned  "${OUT_C4}/results_kept_edges.csv" \
        --scores-cond        "${OUT_C4}/results_scores_cond.npy" \
        --scores-raw         "${OUT_C4}/results_scores_raw.npy" \
        --node-index         "${OUT_C4}/results_node_index.json" \
        --seed-nodes         "${SEEDS_C4}" \
        --chd-pathways       "${REF_C4}" \
        --min-members        8 \
        --max-members        300 \
        --out-csv            "${OUT_C4}/pathway_scores.csv" \
        --out-json           "${OUT_C4}/pathway_metrics.json" \
        --score-variant      degree_norm

    echo "  [C4/${LABEL}] Baseline comparison..."
    python3 "${SCRIPT_DIR}/baseline_enrichment.py" \
        --edges-merged  "${EDGES_MERGED}" \
        --node-index    "${OUT_C4}/results_node_index.json" \
        --scores-raw    "${OUT_C4}/results_scores_raw.npy" \
        --scores-cond   "${OUT_C4}/results_scores_cond.npy" \
        --bifo-scores   "${OUT_C4}/pathway_scores.csv" \
        --chd-pathways  "${REF_C4}" \
        --seed-nodes    "${SEEDS_C4}" \
        --out-csv       "${OUT_C4}/baseline_comparison.csv" \
        --out-json      "${OUT_C4}/baseline_comparison.json" \
        --min-members   8 --max-members 300
}

echo ""
echo "============================================================"
echo "C4: Reactome pathway-split controls"
echo "============================================================"

run_c4 "notch" \
    "${SCRIPT_DIR}/c4_notch_seed_nodes.txt" \
    "${SCRIPT_DIR}/c4_notch_heldout_nodes.txt" \
    "${SCRIPT_DIR}/c4_notch_pathway_reference.txt"

run_c4 "mapk" \
    "${SCRIPT_DIR}/c4_mapk_seed_nodes.txt" \
    "${SCRIPT_DIR}/c4_mapk_heldout_nodes.txt" \
    "${SCRIPT_DIR}/c4_mapk_pathway_reference.txt"

# C4 summary
echo ""
echo "  [C4 summary] Cross-control comparison..."
python3 - <<PYEOF
import json

def show(label, path):
    try:
        m = json.load(open(path))['metrics']
        p10 = m.get('precision_at_10', m.get('top10_precision', 0))
        e10 = m.get('enrichment_at_10', m.get('top10_enrichment_ratio', float('nan')))
        ri  = m.get('rank_improvement_cond_vs_raw', float('nan'))
        mr  = m.get('mean_ref_rank', m.get('mean_rank_chd_cond', float('nan')))
        print(f"  {label:45s}  P@10={p10:.2f}  enrich={e10:.1f}x  mean_rank={mr:.0f}  rank_imp={ri:+.1f}")
    except Exception as ex:
        print(f"  {label:45s}  ERROR: {ex}")

show("CHD curated — BIFO full",          "${OUT_DIR}/pathway_metrics_full.json")
show("C4/Notch — BIFO full",             "${OUT_DIR}/c4_notch/pathway_metrics.json")
show("C4/MAPK  — BIFO full",             "${OUT_DIR}/c4_mapk/pathway_metrics.json")

# For baseline files (json['methods'] not json['metrics']), use best method
def show_best_baseline(label, path):
    try:
        data = json.load(open(path))
        best = max(data['methods'], key=lambda m: m.get('precision_at_10', 0))
        p10  = best.get('precision_at_10', 0)
        ndcg = best.get('ndcg_at_10', 0)
        ap   = best.get('average_precision', 0)
        method = best.get('method', '?')
        print(f"  {label:45s}  best={method}  P@10={p10:.2f}  NDCG@10={ndcg:.3f}  AP={ap:.3f}")
    except Exception as e:
        print(f"  {label:45s}  ERROR: {e}")

show_best_baseline("CHD curated — best baseline", "${OUT_DIR}/baseline_comparison.json")
show_best_baseline("C4/Notch — best baseline",    "${OUT_DIR}/c4_notch/baseline_comparison.json")
show_best_baseline("C4/MAPK  — best baseline",    "${OUT_DIR}/c4_mapk/baseline_comparison.json")
PYEOF

echo ""
echo "============================================================"
echo "All analyses complete."
echo "============================================================"
