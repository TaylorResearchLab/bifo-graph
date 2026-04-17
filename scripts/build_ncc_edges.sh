#!/usr/bin/env bash
# =============================================================================
# build_ncc_edges.sh
# Stage 2.2: Build NCC/cilia custom pathway membership edges.
#
# Usage:
#   bash build_ncc_edges.sh [COHORT] [NCC_DIR]
#   COHORT  defaults to "chd"
#   NCC_DIR defaults to "./ncc_cilia_pathways"
# Inputs:
#   kf_{cohort}_nodes_clean.csv
#   {NCC_DIR}/*.txt  (20 NCC/cilia gene set files)
# Outputs:
#   kf_{cohort}_ncc_membership_edges.csv  (gene -> NCC pathway edges)
#   kf_{cohort}_ncc_pathway_nodes.csv     (20 NCC_CUSTOM concept nodes)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"
COHORT="${1:-chd}"
NCC_DIR="${2:-}"
if [ -z "$NCC_DIR" ]; then
    if [ -d "./ncc_cilia_pathways" ]; then
        NCC_DIR="./ncc_cilia_pathways"
    elif [ -d "$DATA_DIR/ncc_cilia_pathways" ]; then
        NCC_DIR="$DATA_DIR/ncc_cilia_pathways"
    else
        echo "ERROR: ncc_cilia_pathways/ not found in working dir or repo data/"
        echo "Run setup_workspace.sh first, or: ln -s /path/to/repo/data/ncc_cilia_pathways ."
        exit 1
    fi
fi
if [ ! -d "$NCC_DIR" ]; then
    echo "ERROR: NCC pathway directory not found: $NCC_DIR"
    echo "Download ncc_cilia_pathways.zip from the repository and unzip here."
    exit 1
echo "============================================================"
echo "Building NCC membership edges for cohort: $COHORT"
echo "  NCC dir: $NCC_DIR"
python3 "$PIPELINE_DIR/build_ncc_membership_edges.py" \
    --ncc-dir   "$NCC_DIR" \
    --nodes-csv "kf_${COHORT}_nodes_clean.csv" \
    --out-edges "kf_${COHORT}_ncc_membership_edges.csv" \
    --out-nodes "kf_${COHORT}_ncc_pathway_nodes.csv" \
    --verbose
echo ""
echo "NCC build complete."
echo "Next step: bash merge_files.sh $COHORT"
