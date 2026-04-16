#!/usr/bin/env bash
# =============================================================================
# run_kf_nbl_export.sh
# Stage 1: Export KF-NBL graph data from Neo4j using cypher-shell.
# Same structure as run_kf_chd_export.sh — swap cohort prefix only.
# =============================================================================

set -euo pipefail

USER="${1:-neo4j}"
PASS="${2:-neo4j}"
ADDR="${3:-bolt://localhost:7687}"
SEEDS="${4:-maf001}"   # maf001 (default) | maf01 | "" (original 56-gene)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"

CYPHER_SHELL="${CYPHER_SHELL:-cypher-shell}"
if ! command -v "$CYPHER_SHELL" &>/dev/null; then
    echo "ERROR: cypher-shell not found."
    exit 1
fi

run_query() {
    local qfile="$1"
    local outfile="$2"
    echo "[$(date +%H:%M:%S)] Running $qfile -> $outfile ..."
    "$CYPHER_SHELL" -u "$USER" -p "$PASS" -a "$ADDR" \
        --format plain \
        --file "./$qfile" \
        > "$outfile" || {
            if [ ! -s "$outfile" ]; then
                echo "ERROR: $outfile is empty."
                exit 1
            fi
        }
    echo "  Done: $(wc -l < "$outfile") lines written"
}

echo "============================================================"
echo "KF-NBL Neo4j Export"
echo "  Server : $ADDR"
echo "  User   : $USER"
echo "============================================================"


# Generate cypher query files from current seed list
SEEDS_FILE="kf_nbl_seeds${SEEDS:+_${SEEDS}}.txt"
# Resolve seed file: working dir first, then repo data/cohorts/
if [ ! -f "$SEEDS_FILE" ]; then
    REPO_SEEDS="$DATA_DIR/cohorts/nbl/$SEEDS_FILE"
    if [ -f "$REPO_SEEDS" ]; then
        SEEDS_FILE="$REPO_SEEDS"
    else
        echo "ERROR: Seed file not found: $SEEDS_FILE"
        echo "Looked in: $(pwd)/$SEEDS_FILE and $REPO_SEEDS"
        exit 1
    fi
fi
echo "Generating query files from: $SEEDS_FILE"
python3 "$PIPELINE_DIR/generate_export_cypher.py" \
    --seeds   "$SEEDS_FILE" \
    --cohort  "nbl" \
    --out-dir "."
echo ""

run_query kf_nbl_query2.cypher kf_nbl_seed_nodes.csv
run_query kf_nbl_query3.cypher kf_nbl_edges_raw.csv
run_query kf_nbl_query4.cypher kf_nbl_nodes.csv
run_query kf_nbl_query5.cypher kf_nbl_pathway_membership_edges.csv

echo ""
echo "Export complete."
echo "Next step: bash clean_files.sh nbl"
