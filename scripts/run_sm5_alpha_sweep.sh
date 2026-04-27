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
LOGS_DIR="/mnt/isilon/taylor_lab/data/projects/BIFO_2026/logs/sm5_alpha"
mkdir -p "$LOGS_DIR"

cd "$BIFO_GRAPH"

# --- α values to sweep --------------------------------------------------------
ALPHAS=(0.20 0.35 0.50 0.65 0.80)

# --- Pipeline inputs ----------------------------------------------------------
# v0.3.1: Both conditioning and scoring use the cleaned, deduplicated nodes
# file (nodes_clean_noncc.csv.gz). Earlier versions used kf_*_nodes_extended.csv
# for conditioning, which contains duplicate node_id rows that triggered a
# silent dimension mismatch in bifo_conditioning.py outputs. See bifo_conditioning.py
# v0.3.1 patch notes for full diagnosis.
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
for f in "$CHD_SEED_CUIS" "$CHD_REFERENCE" \
         "$CHD_NODES_CLEAN" "$CHD_EDGES_RAW" "$CHD_EDGES_MERGED" \
         "$NBL_SEED_CUIS" "$NBL_REFERENCE" \
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
    local outdir nodes_cond edges_cond seed_cuis reference
    local nodes_clean edges_raw edges_merged
    local disease n_probands n_seeds

    if [[ "$cohort" == "chd" ]]; then
        outdir="$BIFO_GRAPH/results/kf_chd/sensitivity/alpha${alpha}"
        # NOTE v0.3.1: conditioning uses the cleaned, deduplicated nodes file
        # (nodes_clean_noncc.csv.gz), not the raw nodes_extended.csv. The
        # extended file contains duplicate node_id rows that previously caused
        # silent dimension mismatch between scores_cond.npy and node_index.json.
        # bifo_conditioning.py v0.3.1+ also defends against this with a dedup
        # step + assertion, but we pass the clean file here for safety and
        # parity with the canonical scoring pipeline.
        nodes_cond="$CHD_NODES_CLEAN"; edges_cond="$CHD_EDGES_RAW"
        seed_cuis="$CHD_SEED_CUIS"; reference="$CHD_REFERENCE"
        nodes_clean="$CHD_NODES_CLEAN"; edges_raw="$CHD_EDGES_RAW"
        edges_merged="$CHD_EDGES_MERGED"
        disease="congenital_heart_disease"; n_probands=697; n_seeds=1276
    else
        outdir="$BIFO_GRAPH/results/kf_nbl/sensitivity/alpha${alpha}"
        nodes_cond="$NBL_NODES_CLEAN"; edges_cond="$NBL_EDGES_RAW"
        seed_cuis="$NBL_SEED_CUIS"; reference="$NBL_REFERENCE"
        nodes_clean="$NBL_NODES_CLEAN"; edges_raw="$NBL_EDGES_RAW"
        edges_merged="$NBL_EDGES_MERGED"
        disease="neuroblastoma"; n_probands=460; n_seeds=1395
    fi
    mkdir -p "$outdir"

    local log_cond="$LOGS_DIR/${cohort}_alpha${alpha}_conditioning.log"
    local log_score="$LOGS_DIR/${cohort}_alpha${alpha}_scoring.log"
    local log_summ="$LOGS_DIR/${cohort}_alpha${alpha}_summarize.log"

    echo "=== KF-${cohort_upper}  alpha=${alpha}  ==="
    echo "    Outdir: $outdir"

    # ---- 1. bifo_conditioning.py at this α ---
    echo "[1/3] bifo_conditioning.py --alpha $alpha ..."
    python3 "$BIFO_GRAPH/pipeline/bifo_conditioning.py" \
        --nodes         "$nodes_cond" \
        --edges         "$edges_cond" \
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
