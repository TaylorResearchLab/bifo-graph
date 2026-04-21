#!/usr/bin/env bash
# =============================================================================
# run_scoring.sh
# Stage 4: Pathway scoring — standard MSIGDB universe.
#
# Usage:
#   bash run_scoring.sh [COHORT]
#   COHORT defaults to "chd"
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
N_CORES="${2:-0}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
RESULTS_DIR="kf_${COHORT}_results"
DATA_DIR="$REPO_DIR/data"

echo "============================================================"
echo "Pathway Scoring for cohort: $COHORT"
echo "============================================================"

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
    --chd-pathways      "$DATA_DIR/cohorts/${COHORT}/kf_${COHORT}_cilia_reference.txt" \
    --allowed-name-prefixes HALLMARK_ REACTOME_ WP_ KEGG_ BIOCARTA_ PID_ \
    --min-members 8 --max-members 300 \
    --n-permutations 1000 \
    --null-type membership-rewiring \
    --out-csv           "${RESULTS_DIR}/pathway_scores_standard.csv" \
    --out-json          "${RESULTS_DIR}/pathway_metrics_standard.json" \
    --n-cores           "$N_CORES"

echo ""
echo "Scoring complete."
echo "Next step: bash run_baseline.sh $COHORT"
