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
#
# Outputs:
#   kf_{cohort}_nodes_clean.csv
#   kf_{cohort}_edges_raw_clean.csv
#   kf_{cohort}_pathway_membership_edges_clean.csv
#
# What clean_cypher_output.py does:
#   - Strips leading/trailing spaces from column names
#   - Removes double-quotes wrapping field values (cypher-shell artifact)
#   - Skips Java warning lines at the top of the file
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "============================================================"
echo "Cleaning cypher-shell output for cohort: $COHORT"
echo "============================================================"

python3 "$SCRIPT_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_nodes.csv" \
    "kf_${COHORT}_nodes_clean.csv"

python3 "$SCRIPT_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_edges_raw.csv" \
    "kf_${COHORT}_edges_raw_clean.csv"

python3 "$SCRIPT_DIR/clean_cypher_output.py" \
    "kf_${COHORT}_pathway_membership_edges.csv" \
    "kf_${COHORT}_pathway_membership_edges_clean.csv"

echo ""
echo "Cleaning complete."
echo "Next step: bash build_ncc_edges.sh $COHORT"
