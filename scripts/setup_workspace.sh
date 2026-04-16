#!/usr/bin/env bash
# =============================================================================
# setup_workspace.sh
# Prepare a working directory for BIFO pipeline execution.
#
# Usage:
#   bash /path/to/bifo-graph/scripts/setup_workspace.sh [COHORT] [SEEDS] [WORK_DIR]
#
#   COHORT    chd or nbl (default: chd)
#   SEEDS     maf001 (default) | maf01 | "" (original 56-gene list)
#   WORK_DIR  target working directory (default: ./bifo_run_{COHORT})
#
# What it does:
#   Symlinks all required repo data files into WORK_DIR so the pipeline
#   scripts can find them. Run once before running the pipeline.
#
# After running:
#   cd WORK_DIR
#   bash /path/to/bifo-graph/scripts/run_full_pipeline.sh COHORT SEEDS [user] [pass] [addr]
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SEEDS="${2:-maf001}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"
WORK_DIR="${3:-./bifo_run_${COHORT}}"

mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"

echo "============================================================"
echo "Setting up BIFO workspace"
echo "  Repo    : $REPO_DIR"
echo "  Cohort  : $COHORT"
echo "  Seeds   : ${SEEDS:-original (56 genes)}"
echo "  Work dir: $WORK_DIR"
echo "============================================================"
echo ""

link_item() {
    local src="$1"
    local dst="$2"
    if [ ! -e "$src" ]; then
        echo "  WARNING: not found: $src"
        return
    fi
    ln -sfn "$src" "$dst"
    echo "  linked: $(basename "$dst")"
}

# Config
link_item "$REPO_DIR/config/bifo_ddkg_mapping.yaml" "$WORK_DIR/bifo_ddkg_mapping.yaml"

# NCC pathway gene sets
link_item "$REPO_DIR/data/ncc_cilia_pathways" "$WORK_DIR/ncc_cilia_pathways"

# Seed files — link both the versioned and plain names
SEED_SUFFIX="${SEEDS:+_${SEEDS}}"
SEEDS_FILE="$REPO_DIR/data/cohorts/${COHORT}/kf_${COHORT}_seeds${SEED_SUFFIX}.txt"
link_item "$SEEDS_FILE" "$WORK_DIR/kf_${COHORT}_seeds${SEED_SUFFIX}.txt"
if [ -n "$SEEDS" ]; then
    link_item "$SEEDS_FILE" "$WORK_DIR/kf_${COHORT}_seeds.txt"
fi

# NCC reference pathways
link_item "$REPO_DIR/data/cohorts/${COHORT}/kf_${COHORT}_ncc_reference.txt" \
          "$WORK_DIR/kf_${COHORT}_ncc_reference.txt"

echo ""
echo "Workspace ready: $WORK_DIR"
echo ""
echo "Next steps:"
echo "  cd $WORK_DIR"
echo "  bash $SCRIPT_DIR/run_full_pipeline.sh $COHORT ${SEEDS} neo4j PASSWORD bolt://localhost:7687"
