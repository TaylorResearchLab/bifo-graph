#!/usr/bin/env bash
# =============================================================================
# run_kf_chd_export.sh
# Stage 1: Export KF-CHD graph data from Neo4j using cypher-shell.
#
# Usage:
#   bash run_kf_chd_export.sh [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR]
#
# Defaults:
#   NEO4J_USER = neo4j
#   NEO4J_PASS = neo4j
#   NEO4J_ADDR = bolt://localhost:7687
#
# Outputs (written to current directory):
#   kf_chd_seed_nodes.csv
#   kf_chd_nodes.csv
#   kf_chd_edges_raw.csv
#   kf_chd_pathway_membership_edges.csv
#
# Requirements:
#   - cypher-shell on PATH (or set CYPHER_SHELL below)
#   - kf_chd_query2.cypher ... kf_chd_query5.cypher in same directory
#     (generated from kf_chd_export_queries.cypher)
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

# Locate cypher-shell
CYPHER_SHELL="${CYPHER_SHELL:-cypher-shell}"
if ! command -v "$CYPHER_SHELL" &>/dev/null; then
    echo "ERROR: cypher-shell not found. Set CYPHER_SHELL env var or add to PATH."
    exit 1
fi

run_query() {
    local qfile="$1"
    local outfile="$2"
    echo "[$(date +%H:%M:%S)] Running $qfile -> $outfile ..."
    "$CYPHER_SHELL" -u "$USER" -p "$PASS" -a "$ADDR" \
        --format plain \
        < "./$qfile" \
        > "$outfile" 2>/dev/null || {
            # cypher-shell may exit non-zero on warnings — check output has content
            if [ ! -s "$outfile" ]; then
                echo "ERROR: $outfile is empty. Check Neo4j connection and query."
                exit 1
            fi
        }
    local rows
    rows=$(wc -l < "$outfile")
    echo "  Done: $rows lines written"
    if [ "$rows" -le 1 ]; then
        echo "  WARNING: Only $rows line(s) — query may have returned no results."
    fi
}

echo "============================================================"
echo "KF-CHD Neo4j Export"
echo "  Server : $ADDR"
echo "  User   : $USER"
echo "============================================================"


# Generate cypher query files from current seed list
SEEDS_FILE="kf_chd_seeds${SEEDS:+_${SEEDS}}.txt"
if [ ! -f "$SEEDS_FILE" ]; then
    echo "ERROR: Seed file not found: $SEEDS_FILE"
    echo "Expected in working directory. Copy from repo: data/cohorts/chd/$SEEDS_FILE"
    exit 1
fi
echo "Generating query files from: $SEEDS_FILE"
python3 "$PIPELINE_DIR/generate_export_cypher.py" \
    --seeds   "$SEEDS_FILE" \
    --cohort  "chd" \
    --out-dir "."
echo ""

# Query 2: seed concept nodes (small, ~56 rows)
run_query kf_chd_query2.cypher kf_chd_seed_nodes.csv

# Query 3: 1-hop edges
run_query kf_chd_query3.cypher kf_chd_edges_raw.csv

# Query 4: all nodes with HGNC-preferred SAB
run_query kf_chd_query4.cypher kf_chd_nodes.csv

# Query 5: pathway membership edges (largest — ~1.1M rows)
run_query kf_chd_query5.cypher kf_chd_pathway_membership_edges.csv

echo ""
echo "Export complete."
echo "Next step: bash clean_files.sh"
