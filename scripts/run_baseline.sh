#!/usr/bin/env bash
# =============================================================================
# run_baseline.sh
# Stage 5: Baseline enrichment comparison (B1-B4).
#
# Usage:
#   bash run_baseline.sh [COHORT]
#   COHORT defaults to "chd"
#
# Inputs:
#   kf_{cohort}_edges_all.csv
#   kf_{cohort}_results/results_node_index.json
#   kf_{cohort}_results/results_scores_raw.npy
#   kf_{cohort}_results/results_scores_cond.npy
#   kf_{cohort}_results/pathway_scores_standard.csv
#   kf_{cohort}_ncc_reference.txt
#   kf_{cohort}_seed_cuis.txt
#
# Outputs (written to kf_{cohort}_results/):
#   baseline_comparison.csv
#   baseline_comparison.json
#
# Methods run:
#   B1  seed_fisher            Hypergeometric on seed gene set membership
#   B2  neighborhood_fisher    Hypergeometric on 1-hop neighborhood
#   B3  raw_ppr_gsea           Pre-ranked GSEA on raw PPR scores
#   B3b cond_ppr_gsea          Pre-ranked GSEA on conditioned PPR scores
#   B4  bifo_full              BIFO full-arm degree-normalized score
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"
RESULTS_DIR="kf_${COHORT}_results"

echo "============================================================"
echo "Baseline Comparison for cohort: $COHORT"
echo "============================================================"

NCC_REF="kf_${COHORT}_ncc_reference.txt"
if [ ! -f "$NCC_REF" ]; then
    REPO_REF="$DATA_DIR/cohorts/${COHORT}/$NCC_REF"
    if [ -f "$REPO_REF" ]; then
        NCC_REF="$REPO_REF"
    fi
fi

python3 "$PIPELINE_DIR/baseline_enrichment.py" \
    --edges-merged  "kf_${COHORT}_edges_all.csv" \
    --node-index    "${RESULTS_DIR}/results_node_index.json" \
    --scores-raw    "${RESULTS_DIR}/results_scores_raw.npy" \
    --scores-cond   "${RESULTS_DIR}/results_scores_cond.npy" \
    --bifo-scores   "${RESULTS_DIR}/pathway_scores_standard.csv" \
    --chd-pathways  "${NCC_REF:-kf_${COHORT}_ncc_reference.txt}" \
    --seed-nodes    "kf_${COHORT}_seed_cuis.txt" \
    --out-csv       "${RESULTS_DIR}/baseline_comparison.csv" \
    --out-json      "${RESULTS_DIR}/baseline_comparison.json"

echo ""
echo "Baseline comparison complete."
echo "Next step: bash run_resampling.sh $COHORT"
