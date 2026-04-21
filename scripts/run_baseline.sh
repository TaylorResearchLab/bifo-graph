#!/usr/bin/env bash
# =============================================================================
# run_baseline.sh
# Stage 5: Baseline enrichment comparison (B0-B4).
#
# Usage:
#   bash run_baseline.sh [COHORT]
#   COHORT defaults to "chd"
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
RESULTS_DIR="kf_${COHORT}_results"

echo "============================================================"
echo "Baseline Comparison for cohort: $COHORT"
echo "============================================================"

python3 "$PIPELINE_DIR/baseline_enrichment.py" \
    --edges-merged  "kf_${COHORT}_edges_all.csv" \
    --node-index    "${RESULTS_DIR}/results_node_index.json" \
    --scores-raw    "${RESULTS_DIR}/results_scores_raw.npy" \
    --scores-cond   "${RESULTS_DIR}/results_scores_cond.npy" \
    --bifo-scores   "${RESULTS_DIR}/pathway_scores_standard.csv" \
    --chd-pathways  "$DATA_DIR/cohorts/${COHORT}/kf_${COHORT}_cilia_reference.txt" \
    --seed-nodes    "kf_${COHORT}_seed_cuis.txt" \
    --small-universe \
    --out-csv       "${RESULTS_DIR}/baseline_comparison.csv" \
    --out-json      "${RESULTS_DIR}/baseline_comparison.json"

echo ""
echo "Baseline comparison complete."
echo "Next step: bash run_resampling.sh $COHORT"
