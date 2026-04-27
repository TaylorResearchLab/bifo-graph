#!/usr/bin/env bash
# =============================================================================
# run_sm5_alpha_sweep.sh
# §SM5.1 — PPR teleport probability sensitivity analysis
#
# Sweeps α ∈ {0.20, 0.35, 0.50, 0.65, 0.80} for both KF-CHD and KF-NBL.
# Canonical α = 0.50 sits at the center.
#
# Per α value, runs the full pipeline:
#   1. bifo_conditioning.py --alpha α
#   2. score_pathways.py --ppr-alpha α
#   3. summarize_results.py
#
# Pairs CHD + NBL across the 192-core HPC (96 cores each, run concurrently).
# Outputs per α: bifo-graph/results/kf_{cohort}/sensitivity/alpha{α}/
#
# Usage:
#   nohup bash scripts/run_sm5_alpha_sweep.sh \
#     > /mnt/isilon/taylor_lab/data/projects/BIFO_2026/logs/sm5_alpha_master.log 2>&1 &
#
# Determinism:
#   --perm-random-seed 42 across all α points (per pre-specified §SM5 design)
# =============================================================================

set -euo pipefail

# --- Paths --------------------------------------------------------------------
BIFO_GRAPH="/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo-graph"
BIFO_RUN_CHD="/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_chd"
BIFO_RUN_NBL="/mnt/isilon/taylor_lab/data/projects/BIFO_2026/bifo_run_nbl"
LOGS_DIR="/mnt/isilon/taylor_lab/data/projects/BIFO_2026/logs/sm5_alpha"
mkdir -p "$LOGS_DIR"

cd "$BIFO_GRAPH"

# --- α values to sweep --------------------------------------------------------
ALPHAS=(0.20 0.35 0.50 0.65 0.80)

# --- Canonical conditioning inputs (raw graph files) -------------------------
CHD_NODES_EXT="$BIFO_RUN_CHD/kf_chd_nodes_extended.csv"
CHD_EDGES_ALL="$BIFO_RUN_CHD/kf_chd_edges_all.csv"
NBL_NODES_EXT="$BIFO_RUN_NBL/kf_nbl_nodes_extended.csv"
NBL_EDGES_ALL="$BIFO_RUN_NBL/kf_nbl_edges_all.csv"

# --- score_pathways.py inputs (cleaned/conditioned versions) -----------------
CHD_SEED_CUIS="$BIFO_GRAPH/data/cohorts/chd/kf_chd_seed_cuis.txt"
CHD_REFERENCE="$BIFO_GRAPH/data/cohorts/chd/kf_chd_cilia_reference.txt"
CHD_NODES_CLEAN="$BIFO_GRAPH/results/kf_chd/nodes_clean_noncc.csv.gz"
CHD_EDGES_RAW="$BIFO_GRAPH/results/kf_chd/edges_all_noncc.csv.gz"
CHD_EDGES_MERGED="$BIFO_GRAPH/results/kf_chd/edges_all_noncc.csv"

NBL_SEED_CUIS="$BIFO_GRAPH/data/cohorts/nbl/kf_nbl_seed_cuis.txt"
NBL_REFERENCE="$BIFO_GRAPH/data/cohorts/nbl/kf_nbl_cilia_reference.txt"
NBL_NODES_CLEAN="$BIFO_GRAPH/results/kf_nbl/nodes_clean_noncc.csv.gz"
NBL_EDGES_RAW="$BIFO_GRAPH/results/kf_nbl/edges_all_noncc.csv.gz"
NBL_EDGES_MERGED="$BIFO_GRAPH/results/kf_nbl/edges_all_noncc.csv"

MAPPING="$BIFO_GRAPH/config/bifo_mapping.yaml"

# --- Verify all inputs exist before launching anything -----------------------
echo "Verifying input files..."
for f in "$CHD_NODES_EXT" "$CHD_EDGES_ALL" "$CHD_SEED_CUIS" "$CHD_REFERENCE" \
         "$CHD_NODES_CLEAN" "$CHD_EDGES_RAW" "$CHD_EDGES_MERGED" \
         "$NBL_NODES_EXT" "$NBL_EDGES_ALL" "$NBL_SEED_CUIS" "$NBL_REFERENCE" \
         "$NBL_NODES_CLEAN" "$NBL_EDGES_RAW" "$NBL_EDGES_MERGED" \
         "$MAPPING"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: missing input file: $f"
        exit 1
    fi
done
echo "  All input files verified."

# --- Run one cohort at one α --------------------------------------------------
run_cohort_at_alpha() {
    local cohort="$1"
    local alpha="$2"
    local n_cores="$3"

    local cohort_upper="$(echo "$cohort" | tr '[:lower:]' '[:upper:]')"
    local outdir="$BIFO_GRAPH/results/kf_${cohort}/sensitivity/alpha${alpha}"
    mkdir -p "$outdir"

    local nodes_ext edges_all seed_cuis reference nodes_clean edges_raw edges_merged disease n_probands n_seeds
    if [ "$cohort" = "chd" ]; then
        nodes_ext="$CHD_NODES_EXT"; edges_all="$CHD_EDGES_ALL"
        seed_cuis="$CHD_SEED_CUIS"; reference="$CHD_REFERENCE"
        nodes_clean="$CHD_NODES_CLEAN"; edges_raw="$CHD_EDGES_RAW"
        edges_merged="$CHD_EDGES_MERGED"
        disease="congenital heart disease"; n_probands=697; n_seeds=1276
    else
        nodes_ext="$NBL_NODES_EXT"; edges_all="$NBL_EDGES_ALL"
        seed_cuis="$NBL_SEED_CUIS"; reference="$NBL_REFERENCE"
        nodes_clean="$NBL_NODES_CLEAN"; edges_raw="$NBL_EDGES_RAW"
        edges_merged="$NBL_EDGES_MERGED"
        disease="neuroblastoma"; n_probands=460; n_seeds=1395
    fi

    local log_cond="$LOGS_DIR/${cohort}_alpha${alpha}_conditioning.log"
    local log_score="$LOGS_DIR/${cohort}_alpha${alpha}_scoring.log"
    local log_summ="$LOGS_DIR/${cohort}_alpha${alpha}_summarize.log"

    echo "=== KF-${cohort_upper}  alpha=${alpha}  ==="
    echo "    Outdir: $outdir"

    # ---- 1. bifo_conditioning.py at this α ---
    echo "[1/3] bifo_conditioning.py --alpha $alpha ..."
    python3 "$BIFO_GRAPH/pipeline/bifo_conditioning.py" \
        --nodes         "$nodes_ext" \
        --edges         "$edges_all" \
        --mapping       "$MAPPING" \
        --seed-nodes    "$seed_cuis" \
        --heldout-nodes "$seed_cuis" \
        --alpha         "$alpha" \
        --out-json      "$outdir/results.json" \
        > "$log_cond" 2>&1
    echo "      conditioning complete"

    # ---- 2. score_pathways.py at this α ---
    echo "[2/3] score_pathways.py --ppr-alpha $alpha ..."
    python3 "$BIFO_GRAPH/pipeline/score_pathways.py" \
        --nodes             "$nodes_clean" \
        --edges-raw         "$edges_raw" \
        --edges-conditioned "$outdir/results_kept_edges.csv" \
        --scores-cond       "$outdir/results_scores_cond.npy" \
        --scores-raw        "$outdir/results_scores_raw.npy" \
        --node-index        "$outdir/results_node_index.json" \
        --seed-nodes        "$seed_cuis" \
        --chd-pathways      "$reference" \
        --allowed-name-prefixes HALLMARK_ REACTOME_ WP_ KEGG_ BIOCARTA_ PID_ \
        --min-members 8 --max-members 300 \
        --n-permutations 1000 --null-type membership-rewiring \
        --perm-random-seed 42 --n-cores "$n_cores" \
        --ppr-alpha         "$alpha" \
        --out-csv           "$outdir/pathway_scores_standard.csv" \
        --out-json          "$outdir/pathway_metrics_standard.json" \
        --out-member-scores     "$outdir/pathway_member_scores.tsv" \
        --out-influential-nodes "$outdir/pathway_influential_nodes.tsv" \
        > "$log_score" 2>&1
    echo "      scoring complete"

    # ---- 3. summarize_results.py ---
    echo "[3/3] summarize_results.py ..."
    python3 "$BIFO_GRAPH/pipeline/summarize_results.py" \
        --scores            "$outdir/pathway_scores_standard.csv" \
        --seeds             "$seed_cuis" \
        --edges-merged      "$edges_merged" \
        --nodes             "$nodes_clean" \
        --reference         "$reference" \
        --member-scores     "$outdir/pathway_member_scores.tsv" \
        --influential-nodes "$outdir/pathway_influential_nodes.tsv" \
        --cohort-name       "KF-${cohort_upper} alpha=${alpha}" \
        --disease           "$disease" \
        --n-probands        "$n_probands" \
        --n-resolved-seeds  "$n_seeds" \
        --outdir            "$outdir" \
        > "$log_summ" 2>&1
    echo "      summarize complete"

    # ---- 4. Print WP_CILIOPATHIES line for live monitoring ---
    python3 -c "
import csv
for r in csv.DictReader(open('$outdir/pathway_scores_standard.csv')):
    if 'CILIO' in r.get('pathway_name',''):
        print(f\"  WP_CILIOPATHIES @ alpha=$alpha (KF-${cohort_upper}): rank={r['rank']}, null_z={r['null_z']}, q={r['empirical_q']}\")
        break
" 2>&1 || echo "  (WP_CILIOPATHIES not found — investigate $outdir)"
}

# --- Main loop: paired CHD + NBL per α ----------------------------------------
echo ""
echo "============================================================"
echo "§SM5.1 alpha sweep: 5 alpha values x 2 cohorts = 10 full pipeline runs"
echo "Pairs run concurrently (96 cores each, 192 total)"
echo "Started: $(date)"
echo "============================================================"

for alpha in "${ALPHAS[@]}"; do
    echo ""
    echo "######################################################################"
    echo "###  alpha = $alpha  (paired CHD + NBL)  $(date)"
    echo "######################################################################"

    (run_cohort_at_alpha chd "$alpha" 96) > "$LOGS_DIR/pair_chd_alpha${alpha}.log" 2>&1 &
    PID_CHD=$!
    echo "  CHD PID: $PID_CHD"

    (run_cohort_at_alpha nbl "$alpha" 96) > "$LOGS_DIR/pair_nbl_alpha${alpha}.log" 2>&1 &
    PID_NBL=$!
    echo "  NBL PID: $PID_NBL"

    echo "  Waiting on PIDs $PID_CHD and $PID_NBL..."
    wait $PID_CHD; STATUS_CHD=$?
    wait $PID_NBL; STATUS_NBL=$?

    if [ $STATUS_CHD -ne 0 ] || [ $STATUS_NBL -ne 0 ]; then
        echo ""
        echo "ERROR: at least one cohort failed at alpha=$alpha"
        echo "  CHD exit: $STATUS_CHD"
        echo "  NBL exit: $STATUS_NBL"
        echo "  See logs in $LOGS_DIR"
        exit 1
    fi

    grep "WP_CILIOPATHIES" "$LOGS_DIR/pair_chd_alpha${alpha}.log" 2>/dev/null || true
    grep "WP_CILIOPATHIES" "$LOGS_DIR/pair_nbl_alpha${alpha}.log" 2>/dev/null || true
done

echo ""
echo "============================================================"
echo "§SM5.1 alpha sweep COMPLETE  $(date)"
echo "============================================================"
echo ""
echo "Next steps:"
echo "  1. §SM5.3 calibration:  python3 scripts/run_sm5_calibration_sweep.py"
echo "  2. §SM5.2 seed scatter: (separate script, draft tomorrow)"
echo "  3. Aggregate Table SM5.1: (separate script)"
