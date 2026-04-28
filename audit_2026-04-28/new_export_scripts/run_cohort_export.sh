#!/usr/bin/env bash
# =============================================================================
# run_cohort_export.sh
# Stage 1: Export <COHORT> graph data from Neo4j using cypher-shell.
#
# This is the generic per-cohort runner. The five cohort-specific wrappers
# (run_kf_chd_export.sh, run_kf_nbl_export.sh, run_kf_tall_export.sh,
# run_cbtn_export.sh, run_rsbd_export.sh) call this with their cohort name.
#
# Usage:
#   bash run_cohort_export.sh COHORT [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR] [SEEDS]
#
#   COHORT      chd | nbl | tall | cbtn | rsbd
#   NEO4J_USER  default: neo4j
#   NEO4J_PASS  default: neo4j
#   NEO4J_ADDR  default: bolt://localhost:7687
#   SEEDS       default: maf001 (also accepts maf01, maf001_n2, maf001_n3)
#
# Outputs (written to current directory):
#   kf_<cohort>_seed_nodes.csv
#   kf_<cohort>_edges_raw.csv
#   kf_<cohort>_nodes.csv
#   kf_<cohort>_pathway_membership_edges.csv
#   kf_<cohort>_pathway_member_nodes.csv         <- NEW (audit_2026-04-28)
#
# Requirements:
#   - cypher-shell on PATH (or set CYPHER_SHELL env var)
#   - data/cohorts/<cohort>/kf_<cohort>_seeds_<SEEDS>.txt in repo
#   - pipeline/generate_export_cypher.py (patched audit_2026-04-28 version)
#
# CHANGELOG (audit_2026-04-28):
#   - Added Query 6 invocation for kf_<cohort>_pathway_member_nodes.csv.
#     Without this, clean_files.sh and merge_files.sh fail (they expect this
#     file to exist). The previous run_kf_{chd,nbl}_export.sh scripts only ran
#     queries 2-5 and the missing Q6 output file had to be supplied manually.
#   - Generalized for five cohorts (chd, nbl, tall, cbtn, rsbd).
# =============================================================================

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "ERROR: COHORT argument required."
    echo "Usage: bash run_cohort_export.sh COHORT [NEO4J_USER] [NEO4J_PASS] [NEO4J_ADDR] [SEEDS]"
    exit 1
fi

COHORT="$1"
USER="${2:-neo4j}"
PASS="${3:-neo4j}"
ADDR="${4:-bolt://localhost:7687}"
SEEDS="${5:-maf001}"

# Validate cohort
case "$COHORT" in
    chd|nbl|tall|cbtn|rsbd) ;;
    *) echo "ERROR: Unknown cohort '$COHORT'. Expected one of: chd, nbl, tall, cbtn, rsbd"
       exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"

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
        --file "./$qfile" \
        > "$outfile" || {
            # cypher-shell may exit non-zero on warnings -- check output has content
            if [ ! -s "$outfile" ]; then
                echo "ERROR: $outfile is empty. Check Neo4j connection and query."
                exit 1
            fi
        }
    local rows
    rows=$(wc -l < "$outfile")
    echo "  Done: $rows lines written"
    if [ "$rows" -le 1 ]; then
        echo "  WARNING: Only $rows line(s) -- query may have returned no results."
    fi
}

echo "============================================================"
echo "$(echo "$COHORT" | tr '[:lower:]' '[:upper:]') Neo4j Export"
echo "  Server  : $ADDR"
echo "  User    : $USER"
echo "  Seeds   : $SEEDS"
echo "============================================================"

# Resolve seed file: working dir first, then repo data/cohorts/<cohort>/
SEEDS_FILE="kf_${COHORT}_seeds${SEEDS:+_${SEEDS}}.txt"
if [ ! -f "$SEEDS_FILE" ]; then
    REPO_SEEDS="$DATA_DIR/cohorts/${COHORT}/$SEEDS_FILE"
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
    --cohort  "$COHORT" \
    --out-dir "."
echo ""

# Capture md5 of seed file for provenance
SEEDS_MD5=$(md5sum "$SEEDS_FILE" | awk '{print $1}')
echo "Seed file md5: $SEEDS_MD5"
echo ""

# Query 2: seed concept nodes
run_query "kf_${COHORT}_query2.cypher" "kf_${COHORT}_seed_nodes.csv"

# Query 3: 1-hop edges
run_query "kf_${COHORT}_query3.cypher" "kf_${COHORT}_edges_raw.csv"

# Query 4: all nodes (seeds + 1-hop) with HGNC-preferred SAB
run_query "kf_${COHORT}_query4.cypher" "kf_${COHORT}_nodes.csv"

# Query 5: pathway membership edges (largest -- ~1M+ rows typically)
run_query "kf_${COHORT}_query5.cypher" "kf_${COHORT}_pathway_membership_edges.csv"

# Query 6: pathway member nodes (full MSIGDB membership, cohort-independent)
# This is the file added in audit_2026-04-28; without it clean_files.sh fails.
run_query "kf_${COHORT}_query6.cypher" "kf_${COHORT}_pathway_member_nodes.csv"

# Write a small provenance file for this export
PROVENANCE_FILE="kf_${COHORT}_export_provenance.txt"
{
    echo "BIFO Cypher Export Provenance"
    echo "============================="
    echo "Cohort:         $COHORT"
    echo "Date:           $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "Hostname:       $(hostname)"
    echo "Seeds file:     $SEEDS_FILE"
    echo "Seeds md5:      $SEEDS_MD5"
    echo "Seeds count:    $(grep -vc '^#' "$SEEDS_FILE" 2>/dev/null || echo 'unknown')"
    echo "Neo4j addr:     $ADDR"
    echo "Neo4j user:     $USER"
    echo ""
    echo "Output md5s:"
    md5sum "kf_${COHORT}_seed_nodes.csv" \
           "kf_${COHORT}_edges_raw.csv" \
           "kf_${COHORT}_nodes.csv" \
           "kf_${COHORT}_pathway_membership_edges.csv" \
           "kf_${COHORT}_pathway_member_nodes.csv" 2>/dev/null
} > "$PROVENANCE_FILE"

echo ""
echo "Export complete."
echo "Provenance: $PROVENANCE_FILE"
echo "Next step: bash clean_files.sh $COHORT"
