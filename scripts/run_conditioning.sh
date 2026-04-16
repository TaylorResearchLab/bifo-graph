#!/usr/bin/env bash
# =============================================================================
# run_conditioning.sh
# Stage 3: BIFO conditioning and PPR propagation.
#
# Usage:
#   bash run_conditioning.sh [COHORT] [MAPPING]
#   COHORT  defaults to "chd"
#   MAPPING defaults to "./bifo_ddkg_mapping.yaml"
#
# Inputs:
#   kf_{cohort}_nodes_extended.csv
#   kf_{cohort}_edges_all.csv
#   {MAPPING}
#   kf_{cohort}_seed_cuis.txt
#
# Outputs (written to kf_{cohort}_results/):
#   results.json
#   results_kept_edges.csv
#   results_node_index.json
#   results_scores_cond.npy
#   results_scores_raw.npy
#
# Frozen parameters (from CHD curated benchmark):
#   alpha=0.50  (set inside bifo_ddkg_mapping.yaml)
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
MAPPING="${2:-./bifo_ddkg_mapping.yaml}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR="kf_${COHORT}_results"

if [ ! -f "$MAPPING" ]; then
    echo "ERROR: Mapping file not found: $MAPPING"
    exit 1
fi

mkdir -p "$RESULTS_DIR"

echo "============================================================"
echo "BIFO Conditioning for cohort: $COHORT"
echo "  Mapping : $MAPPING"
echo "  Results : $RESULTS_DIR/"
echo "============================================================"

python3 "$SCRIPT_DIR/bifo_conditioning.py" \
    --nodes         "kf_${COHORT}_nodes_extended.csv" \
    --edges         "kf_${COHORT}_edges_all.csv" \
    --mapping       "$MAPPING" \
    --seed-nodes    "kf_${COHORT}_seed_cuis.txt" \
    --heldout-nodes "kf_${COHORT}_seed_cuis.txt" \
    --out-json      "${RESULTS_DIR}/results.json"

echo ""
echo "Conditioning complete."
echo "Next step: bash run_scoring.sh $COHORT"
