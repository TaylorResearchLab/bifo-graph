#!/usr/bin/env bash
# =============================================================================
# run_resampling.sh
# Stage 6: Bootstrap resampling for pathway recovery stability.
#
# Usage:
#   bash run_resampling.sh [COHORT] [N_BOOTS] [N_CORES]
#   COHORT   defaults to "chd"
#   N_BOOTS  defaults to 500 (draws per seed size)
#   N_CORES  defaults to all available cores
#
# Inputs:
#   kf_{cohort}_results/results_kept_edges.csv
#   kf_{cohort}_edges_all.csv
#   kf_{cohort}_results/results_node_index.json
#   kf_{cohort}_results/pathway_scores_standard.csv
#   kf_{cohort}_seeds.txt
#   kf_{cohort}_ncc_reference.txt
#
# Outputs (written to kf_{cohort}_results/):
#   resampling_results.csv    (one row per bootstrap run)
#   resampling_summary.json   (aggregated metrics by seed size)
#
# Design: 500 bootstrap draws at each of 3 seed sizes (10, 20, 30)
#   = 1,500 total runs per cohort
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
N_BOOTS="${2:-500}"
N_CORES="${3:-}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR="kf_${COHORT}_results"

# Default cores to nproc if not specified
if [ -z "$N_CORES" ]; then
    N_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
fi

echo "============================================================"
echo "Bootstrap Resampling for cohort: $COHORT"
echo "  Bootstraps per seed size : $N_BOOTS"
echo "  Seed sizes               : 10, 20, 30"
echo "  Total runs               : $((N_BOOTS * 3))"
echo "  Cores                    : $N_CORES"
echo "============================================================"

python3 "$SCRIPT_DIR/kf_resampling.py" \
    --kept-edges   "${RESULTS_DIR}/results_kept_edges.csv" \
    --edges-merged "kf_${COHORT}_edges_all.csv" \
    --node-index   "${RESULTS_DIR}/results_node_index.json" \
    --bifo-scores  "${RESULTS_DIR}/pathway_scores_standard.csv" \
    --seed-pool    "kf_${COHORT}_seeds.txt" \
    --ref-pathways "kf_${COHORT}_ncc_reference.txt" \
    --out-csv      "${RESULTS_DIR}/resampling_results.csv" \
    --out-json     "${RESULTS_DIR}/resampling_summary.json" \
    --seed-sizes   10 20 30 \
    --n-boots      "$N_BOOTS" \
    --n-cores      "$N_CORES"

echo ""
echo "Resampling complete."
echo "Results: ${RESULTS_DIR}/resampling_summary.json"
