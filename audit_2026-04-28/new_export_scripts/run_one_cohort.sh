#!/usr/bin/env bash
# =============================================================================
# run_one_cohort.sh
# Full BIFO-PPR pipeline for one cohort: Cypher export -> conditioning -> scoring.
#
# Usage:
#   bash run_one_cohort.sh COHORT [SEEDS] [N_CORES] [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR]
#
#   COHORT       chd | nbl | tall | cbtn | rsbd
#   SEEDS        default: maf001
#   N_CORES      default: 0 (auto-detect)
#   NEO4J_USER   default: neo4j
#   NEO4J_PASS   default: neo4j
#   NEO4J_ADDR   default: bolt://localhost:7687
#
# This is the rerun driver. It runs all stages for a single cohort in a fresh
# working directory: _rerun_2026-04-28/<COHORT>/. Logs are tee'd to that dir.
#
# Designed to run in parallel across HPC nodes -- each node calls this with a
# different COHORT argument. All outputs land in the shared filesystem under
# _rerun_2026-04-28/.
# =============================================================================

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: bash run_one_cohort.sh COHORT [SEEDS] [N_CORES] [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR]"
    exit 1
fi

COHORT="$1"
SEEDS="${2:-maf001}"
N_CORES="${3:-0}"
NEO4J_USER="${4:-neo4j}"
NEO4J_PASS="${5:-neo4j}"
NEO4J_ADDR="${6:-bolt://localhost:7687}"

# Validate cohort
case "$COHORT" in
    chd|nbl|tall|cbtn|rsbd) ;;
    *) echo "ERROR: Unknown cohort '$COHORT'. Expected one of: chd, nbl, tall, cbtn, rsbd"
       exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
SCRIPTS_DIR="$REPO_DIR/scripts"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"

# Working directory for this cohort's rerun
RERUN_DIR="$REPO_DIR/_rerun_2026-04-28/$COHORT"
mkdir -p "$RERUN_DIR"

# Log file
LOG_FILE="$RERUN_DIR/run_log.txt"

# Sanity: confirm conda environment is bifo-spectral (Python 3.11)
if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
    echo "Conda environment: $CONDA_DEFAULT_ENV" | tee -a "$LOG_FILE"
    if [ "$CONDA_DEFAULT_ENV" != "bifo-spectral" ]; then
        echo "  WARNING: Expected 'bifo-spectral', got '$CONDA_DEFAULT_ENV'." | tee -a "$LOG_FILE"
        echo "  Run: conda activate bifo-spectral" | tee -a "$LOG_FILE"
    fi
fi
PYTHON_VERSION=$(python3 --version 2>&1)
echo "Python: $PYTHON_VERSION" | tee -a "$LOG_FILE"
if [[ "$PYTHON_VERSION" != *"Python 3.11"* ]]; then
    echo "  WARNING: Expected Python 3.11, got $PYTHON_VERSION" | tee -a "$LOG_FILE"
fi

log_stage() {
    {
        echo ""
        echo "============================================================"
        echo "  STAGE $1: $2  (cohort: $COHORT)"
        echo "  $(date '+%Y-%m-%d %H:%M:%S')"
        echo "============================================================"
    } | tee -a "$LOG_FILE"
}

cd "$RERUN_DIR"

log_stage "1" "Neo4j Cypher export"
bash "$SCRIPT_DIR/run_kf_${COHORT}_export.sh" \
    "$NEO4J_USER" "$NEO4J_PASS" "$NEO4J_ADDR" "$SEEDS" 2>&1 | tee -a "$LOG_FILE"

log_stage "2.1" "Clean cypher-shell output"
bash "$SCRIPTS_DIR/clean_files.sh" "$COHORT" 2>&1 | tee -a "$LOG_FILE"

log_stage "2.2" "Merge CSV files"
bash "$SCRIPTS_DIR/merge_files.sh" "$COHORT" 2>&1 | tee -a "$LOG_FILE"

log_stage "2.3" "Seed CUI lookup"
bash "$SCRIPTS_DIR/run_seed_lookup.sh" "$COHORT" "$SEEDS" 2>&1 | tee -a "$LOG_FILE"

# Capture md5s of the conditioning inputs before running (for provenance)
{
    echo ""
    echo "Conditioning input md5s:"
    md5sum "kf_${COHORT}_nodes_extended.csv" \
           "kf_${COHORT}_edges_all.csv" \
           "kf_${COHORT}_seed_cuis.txt" 2>/dev/null
} | tee -a "$LOG_FILE"

log_stage "3" "BIFO conditioning + PPR"
bash "$SCRIPTS_DIR/run_conditioning.sh" "$COHORT" 2>&1 | tee -a "$LOG_FILE"

log_stage "4" "Pathway scoring"
bash "$SCRIPTS_DIR/run_scoring.sh" "$COHORT" "$N_CORES" 2>&1 | tee -a "$LOG_FILE"

# Capture WP_CILIOPATHIES result for fast comparison
{
    echo ""
    echo "============================================================"
    echo "  WP_CILIOPATHIES result for $COHORT:"
    echo "============================================================"
    SCORES_CSV="kf_${COHORT}_results/pathway_scores_standard.csv"
    if [ -f "$SCORES_CSV" ]; then
        head -1 "$SCORES_CSV"
        grep "WP_CILIOPATHIES" "$SCORES_CSV" || echo "  (WP_CILIOPATHIES not found in scores)"
    else
        echo "  ERROR: Scores file not found at $SCORES_CSV"
    fi
    echo ""
    echo "Pipeline complete for $COHORT."
    echo "Next: optional baseline (run_baseline.sh) and resampling (run_resampling.sh)"
    echo "Outputs in: $RERUN_DIR/kf_${COHORT}_results/"
} | tee -a "$LOG_FILE"
