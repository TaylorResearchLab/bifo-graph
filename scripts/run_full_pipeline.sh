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
# Requirements:
#   Python scripts  : bifo_conditioning.py, score_pathways.py,
#                     baseline_enrichment.py, seed_cui_lookup.py,
#                     clean_cypher_output.py, kf_resampling.py,
#                     summarize_results.py
#   Shell scripts   : run_kf_{cohort}_export.sh, clean_files.sh,
#                     merge_files.sh, run_seed_lookup.sh,
#                     run_conditioning.sh, run_scoring.sh, run_baseline.sh,
#                     run_resampling.sh
#   Data files      : bifo_mapping.yaml, kf_{cohort}_seeds.txt
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SEEDS="${2:-maf001}"
N_CORES="${3:-0}"
NEO4J_USER="${4:-neo4j}"
NEO4J_PASS="${5:-neo4j}"
NEO4J_ADDR="${6:-bolt://localhost:7687}"
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

log_stage "2.2" "Merge CSV files"
bash "$SCRIPT_DIR/merge_files.sh" "$COHORT"

log_stage "2.3" "Seed CUI lookup"
bash "$SCRIPT_DIR/run_seed_lookup.sh" "$COHORT" "$SEEDS"

log_stage "3" "BIFO conditioning + PPR"
bash "$SCRIPT_DIR/run_conditioning.sh" "$COHORT"

log_stage "4" "Pathway scoring"
bash "$SCRIPT_DIR/run_scoring.sh" "$COHORT" "$N_CORES"

log_stage "5" "Baseline enrichment comparison"
bash "$SCRIPT_DIR/run_baseline.sh" "$COHORT"

log_stage "6" "Bootstrap resampling"
bash "$SCRIPT_DIR/run_resampling.sh" "$COHORT" "$SEEDS"

log_stage "7" "Generate summary output files"
RESULTS_DIR="$REPO_DIR/results/kf_${COHORT}"
COHORT_DATA="$DATA_DIR/cohorts/${COHORT}"
COHORT_LABEL="KF-$(echo "$COHORT" | tr '[:lower:]' '[:upper:]')"

# Set disease label per cohort
if [ "$COHORT" = "chd" ]; then
    DISEASE="congenital heart disease"
    N_PROBANDS=697
elif [ "$COHORT" = "nbl" ]; then
    DISEASE="neuroblastoma"
    N_PROBANDS=460
else
    DISEASE=""
    N_PROBANDS=0
fi

python "$PIPELINE_DIR/summarize_results.py" \
    --scores      "$RESULTS_DIR/pathway_scores_standard.csv" \
    --seeds       "$COHORT_DATA/kf_${COHORT}_seeds_maf001.txt" \
    --reference   "$COHORT_DATA/kf_${COHORT}_cilia_reference.txt" \
    --cohort-name "$COHORT_LABEL" \
    --disease     "$DISEASE" \
    --n-probands  "$N_PROBANDS" \
    --outdir      "$RESULTS_DIR"

echo ""
echo "============================================================"
echo "  PIPELINE COMPLETE — cohort: $COHORT"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "  Key outputs:"
echo "    results/kf_${COHORT}/pathway_metrics_standard.json"
echo "    results/kf_${COHORT}/baseline_comparison.json"
echo "    results/kf_${COHORT}/resampling_summary.json"
echo "    results/kf_${COHORT}/pathway_results_summary.tsv"
echo "    results/kf_${COHORT}/pathway_results_llm.md"
echo "============================================================"
