#!/usr/bin/env bash
# =============================================================================
# run_all_cohorts.sh
# Runs the full BIFO-PPR pipeline for all five cohorts.
#
# Usage:
#   bash run_all_cohorts.sh [--parallel] [--cohorts "chd nbl tall cbtn rsbd"] \
#                           [SEEDS] [N_CORES] [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR]
#
# Default behavior:
#   - Runs cohorts sequentially in the order: chd nbl tall cbtn rsbd
#   - Uses SEEDS=maf001, N_CORES=0 (auto), default Neo4j credentials
#
# Examples:
#   # Sequential, all defaults:
#   bash run_all_cohorts.sh
#
#   # Parallel using GNU parallel (one cohort per HPC node, requires shared FS):
#   bash run_all_cohorts.sh --parallel
#
#   # Just CHD and NBL (e.g., for Phase 1 verification):
#   bash run_all_cohorts.sh --cohorts "chd nbl"
#
# Phase 1 recommendation: run CHD and NBL first to verify against canonical,
# THEN run the other three cohorts.
#   bash run_all_cohorts.sh --cohorts "chd nbl"
#   # ... verify results ...
#   bash run_all_cohorts.sh --cohorts "tall cbtn rsbd"
# =============================================================================

set -euo pipefail

PARALLEL=0
COHORTS="chd nbl tall cbtn rsbd"

# Parse flags
while [ $# -gt 0 ]; do
    case "$1" in
        --parallel)
            PARALLEL=1
            shift
            ;;
        --cohorts)
            COHORTS="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

SEEDS="${1:-maf001}"
N_CORES="${2:-0}"
NEO4J_USER="${3:-neo4j}"
NEO4J_PASS="${4:-neo4j}"
NEO4J_ADDR="${5:-bolt://localhost:7687}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
LOG_DIR="$REPO_DIR/_rerun_2026-04-28/_logs"
mkdir -p "$LOG_DIR"

START_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo "============================================================"
echo "BIFO-PPR Five-Cohort Rerun"
echo "  Started:    $START_TIME"
echo "  Cohorts:    $COHORTS"
echo "  Mode:       $([ $PARALLEL -eq 1 ] && echo 'parallel' || echo 'sequential')"
echo "  Seeds:      $SEEDS"
echo "  Cores:      $N_CORES"
echo "  Repo:       $REPO_DIR"
echo "  Logs:       $LOG_DIR"
echo "============================================================"

run_cohort() {
    local cohort="$1"
    local log="$LOG_DIR/${cohort}_$(date +%Y%m%d_%H%M%S).log"
    echo "[$(date +%H:%M:%S)] Starting $cohort ... (log: $log)"
    bash "$SCRIPT_DIR/run_one_cohort.sh" "$cohort" "$SEEDS" "$N_CORES" \
         "$NEO4J_USER" "$NEO4J_PASS" "$NEO4J_ADDR" \
         > "$log" 2>&1 \
         && echo "[$(date +%H:%M:%S)] DONE $cohort" \
         || echo "[$(date +%H:%M:%S)] FAILED $cohort (see $log)"
}

if [ $PARALLEL -eq 1 ]; then
    # Parallel mode: launch all cohorts in background
    pids=()
    for cohort in $COHORTS; do
        run_cohort "$cohort" &
        pids+=($!)
    done

    # Wait for all
    failed=0
    for pid in "${pids[@]}"; do
        wait "$pid" || failed=1
    done

    if [ $failed -ne 0 ]; then
        echo "ERROR: One or more cohort runs failed. Check $LOG_DIR/."
        exit 1
    fi
else
    # Sequential mode
    for cohort in $COHORTS; do
        run_cohort "$cohort"
    done
fi

END_TIME=$(date '+%Y-%m-%d %H:%M:%S')
echo "============================================================"
echo "  Started:  $START_TIME"
echo "  Finished: $END_TIME"
echo "============================================================"

# Summary table of WP_CILIOPATHIES results across cohorts
echo ""
echo "============================================================"
echo "  WP_CILIOPATHIES SUMMARY ACROSS COHORTS"
echo "============================================================"
printf "%-8s  %-12s  %-12s  %-12s\n" "COHORT" "RANK" "NULL_Z" "EMPIRICAL_Q"
echo "----------------------------------------------------------"
for cohort in $COHORTS; do
    SCORES_CSV="$REPO_DIR/_rerun_2026-04-28/${cohort}/kf_${cohort}_results/pathway_scores_standard.csv"
    if [ -f "$SCORES_CSV" ]; then
        # Extract rank, null_z, empirical_q from the row matching WP_CILIOPATHIES
        # Note: rank is computed by sorting on degree_norm (descending), so we
        # compute it inline.
        python3 - <<PYEOF
import pandas as pd
import sys
try:
    df = pd.read_csv("$SCORES_CSV")
    # Sort by degree_norm descending to compute rank
    df_sorted = df.sort_values('degree_norm', ascending=False).reset_index(drop=True)
    df_sorted['rank'] = df_sorted.index + 1
    cilio = df_sorted[df_sorted['name'] == 'WP_CILIOPATHIES']
    if len(cilio) > 0:
        row = cilio.iloc[0]
        print(f"{'$cohort':<8}  {row['rank']:<12}  {row.get('null_z', 'N/A'):<12.3f}  {row.get('empirical_q', 'N/A'):<12.4f}")
    else:
        print(f"{'$cohort':<8}  {'NOT_FOUND':<12}")
except Exception as e:
    print(f"{'$cohort':<8}  ERROR: {e}", file=sys.stderr)
PYEOF
    else
        printf "%-8s  %s\n" "$cohort" "(scores file missing)"
    fi
done

echo "============================================================"
