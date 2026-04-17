#!/usr/bin/env bash
# =============================================================================
# run_scoring.sh
# Stage 4: Pathway scoring — standard universe and NCC/cilia reference.
#
# Usage:
#   bash run_scoring.sh [COHORT]
#   COHORT defaults to "chd"
#
# Inputs:
#   kf_{cohort}_nodes_extended.csv
#   kf_{cohort}_edges_all.csv
#   kf_{cohort}_results/results_kept_edges.csv
#   kf_{cohort}_results/results_scores_cond.npy
#   kf_{cohort}_results/results_scores_raw.npy
#   kf_{cohort}_results/results_node_index.json
#   kf_{cohort}_seed_cuis.txt
#   kf_{cohort}_ncc_reference.txt
#
# Outputs (written to kf_{cohort}_results/):
#   pathway_scores_standard.csv
#   pathway_metrics_standard.json
#   pathway_scores_ncc.csv
#   pathway_metrics_ncc.json
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
N_CORES="${2:-0}"   # 0 = auto-detect all available cores
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"
RESULTS_DIR="kf_${COHORT}_results"

echo "============================================================"
echo "Pathway Scoring for cohort: $COHORT"
echo "============================================================"

# --- Run 4.1: Standard pathway universe (no reference) ---
echo ""
NCC_REF="kf_${COHORT}_ncc_reference.txt"
if [ ! -f "$NCC_REF" ]; then
    REPO_REF="$DATA_DIR/cohorts/${COHORT}/$NCC_REF"
    [ -f "$REPO_REF" ] && NCC_REF="$REPO_REF"
fi

SEED_CUIS="kf_${COHORT}_seed_cuis.txt"

echo "[4.1] Standard universe scoring ..."
python3 "$PIPELINE_DIR/score_pathways.py" \
    --nodes             "kf_${COHORT}_nodes_extended.csv" \
    --edges-raw         "kf_${COHORT}_edges_all.csv" \
    --edges-conditioned "${RESULTS_DIR}/results_kept_edges.csv" \
    --scores-cond       "${RESULTS_DIR}/results_scores_cond.npy" \
    --scores-raw        "${RESULTS_DIR}/results_scores_raw.npy" \
    --node-index        "${RESULTS_DIR}/results_node_index.json" \
    --seed-nodes        "$SEED_CUIS" \
    --out-csv           "${RESULTS_DIR}/pathway_scores_standard.csv" \
    --out-json          "${RESULTS_DIR}/pathway_metrics_standard.json" \
    --n-cores           "$N_CORES"

# --- Run 4.2: NCC/cilia reference evaluation ---
echo ""
echo "[4.2] NCC/cilia reference scoring ..."
python3 "$PIPELINE_DIR/score_pathways.py" \
    --nodes             "kf_${COHORT}_nodes_extended.csv" \
    --edges-raw         "kf_${COHORT}_edges_all.csv" \
    --edges-conditioned "${RESULTS_DIR}/results_kept_edges.csv" \
    --scores-cond       "${RESULTS_DIR}/results_scores_cond.npy" \
    --scores-raw        "${RESULTS_DIR}/results_scores_raw.npy" \
    --node-index        "${RESULTS_DIR}/results_node_index.json" \
    --seed-nodes        "$SEED_CUIS" \
    --chd-pathways      "$NCC_REF" \
    --out-csv           "${RESULTS_DIR}/pathway_scores_ncc.csv" \
    --out-json          "${RESULTS_DIR}/pathway_metrics_ncc.json" \
    --n-cores           "$N_CORES"

echo ""
echo "Scoring complete."
echo "Next step: bash run_baseline.sh $COHORT"
