#!/usr/bin/env bash
# =============================================================================
# run_kf_nbl_export.sh
# Stage 1 wrapper for NBL cohort. Delegates to run_cohort_export.sh.
#
# Usage:
#   bash run_kf_nbl_export.sh [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR] [SEEDS]
#
# Defaults:
#   NEO4J_USER = neo4j
#   NEO4J_PASS = neo4j
#   NEO4J_ADDR = bolt://localhost:7687
#   SEEDS      = maf001
#
# Outputs (written to current directory):
#   kf_nbl_seed_nodes.csv
#   kf_nbl_edges_raw.csv
#   kf_nbl_nodes.csv
#   kf_nbl_pathway_membership_edges.csv
#   kf_nbl_pathway_member_nodes.csv
#   kf_nbl_export_provenance.txt
# =============================================================================

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec bash "$SCRIPT_DIR/run_cohort_export.sh" nbl "$@"
