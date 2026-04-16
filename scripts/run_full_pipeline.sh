#!/usr/bin/env bash
# =============================================================================
# run_full_pipeline.sh
# Master script: run all BIFO pipeline stages end-to-end for one cohort.
#
# Usage:
#   bash run_full_pipeline.sh [COHORT] [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR]
#
#   COHORT      chd (default) or nbl
#   NEO4J_USER  neo4j (default)
#   NEO4J_PASS  neo4j (default)
#   NEO4J_ADDR  bolt://localhost:7687 (default)
#
# Example:
#   bash run_full_pipeline.sh chd neo4j mypassword bolt://localhost:7687
#   bash run_full_pipeline.sh nbl neo4j mypassword bolt://localhost:7687
#
# Requirements (all in same directory as this script):
#   Python scripts  : bifo_conditioning.py, score_pathways.py,
#                     baseline_enrichment.py, seed_cui_lookup.py,
#                     build_ncc_membership_edges.py, clean_cypher_output.py,
#                     kf_resampling.py
#   Shell scripts   : run_kf_{cohort}_export.sh, clean_files.sh,
#                     build_ncc_edges.sh, merge_files.sh, run_seed_lookup.sh,
#                     run_conditioning.sh, run_scoring.sh, run_baseline.sh,
#                     run_resampling.sh
#   Data files      : ncc_cilia_pathways/ directory, bifo_ddkg_mapping.yaml,
#                     kf_{cohort}_seeds.txt, kf_{cohort}_ncc_reference.txt,
#                     kf_{cohort}_query[2-5].cypher
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SEEDS="${2:-maf001}"   # maf001 (default) | maf01 | "" (original)
NEO4J_USER="${3:-neo4j}"
NEO4J_PASS="${4:-neo4j}"
NEO4J_ADDR="${5:-bolt://localhost:7687}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"

log_stage() {
    echo ""
    echo "============================================================"
    echo "  STAGE $1: $2"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "============================================================"
}

log_stage "1" "Neo4j graph export"
bash "$SCRIPT_DIR/run_kf_${COHORT}_export.sh" \
    "$NEO4J_USER" "$NEO4J_PASS" "$NEO4J_ADDR"

log_stage "2.1" "Clean cypher-shell output"
bash "$SCRIPT_DIR/clean_files.sh" "$COHORT"

log_stage "2.2" "Build NCC/cilia membership edges"
bash "$SCRIPT_DIR/build_ncc_edges.sh" "$COHORT"

log_stage "2.3" "Merge CSV files"
bash "$SCRIPT_DIR/merge_files.sh" "$COHORT"

log_stage "2.4" "Seed CUI lookup"
bash "$SCRIPT_DIR/run_seed_lookup.sh" "$COHORT" "$SEEDS"

log_stage "3" "BIFO conditioning + PPR"
bash "$SCRIPT_DIR/run_conditioning.sh" "$COHORT"

log_stage "4" "Pathway scoring (standard + NCC)"
bash "$SCRIPT_DIR/run_scoring.sh" "$COHORT"

log_stage "5" "Baseline enrichment comparison"
bash "$SCRIPT_DIR/run_baseline.sh" "$COHORT"

log_stage "6" "Bootstrap resampling"
bash "$SCRIPT_DIR/run_resampling.sh" "$COHORT" "$SEEDS"

echo ""
echo "============================================================"
echo "  PIPELINE COMPLETE — cohort: $COHORT"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "  Key outputs:"
echo "    kf_${COHORT}_results/pathway_metrics_standard.json"
echo "    kf_${COHORT}_results/pathway_metrics_ncc.json"
echo "    kf_${COHORT}_results/baseline_comparison.json"
echo "    kf_${COHORT}_results/resampling_summary.json"
echo "============================================================"
