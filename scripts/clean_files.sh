#!/usr/bin/env bash
# =============================================================================
# clean_files.sh
# Stage 2.1: Clean raw cypher-shell CSV output files.
#
# Usage:
#   bash clean_files.sh [COHORT]
#   COHORT defaults to "chd"; use "nbl" for neuroblastoma cohort.
#
# Inputs:
#   kf_{cohort}_nodes.csv
#   kf_{cohort}_edges_raw.csv
#   kf_{cohort}_pathway_membership_edges.csv
#   kf_{cohort}_pathway_member_nodes.csv
#
# Outputs:
#   kf_{cohort}_nodes_clean.csv
#   kf_{cohort}_edges_raw_clean.csv
#   kf_{cohort}_pathway_membership_edges_clean.csv
#   kf_{cohort}_pathway_member_nodes_clean.csv
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"

echo "============================================================"
echo "Cleaning cypher-shell output for cohort: $COHORT"
echo "============================================================"

python3 "$PIPELINE_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_nodes.csv" \
    "kf_${COHORT}_nodes_clean.csv"

python3 "$PIPELINE_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_edges_raw.csv" \
    "kf_${COHORT}_edges_raw_clean.csv"

python3 "$PIPELINE_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_pathway_membership_edges.csv" \
    "kf_${COHORT}_pathway_membership_edges_clean.csv"

python3 "$PIPELINE_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_pathway_member_nodes.csv" \
    "kf_${COHORT}_pathway_member_nodes_clean.csv"

echo ""
echo "Cleaning complete."
echo "Next step: bash merge_files.sh $COHORT"
