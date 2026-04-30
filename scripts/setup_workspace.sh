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
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SEEDS="${2:-maf001}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$REPO_DIR/data"
WORK_DIR="${3:-./bifo_run_${COHORT}}"

mkdir -p "$WORK_DIR"
WORK_DIR="$(cd "$WORK_DIR" && pwd)"

echo "============================================================"
echo "Setting up BIFO workspace"
echo "  Repo    : $REPO_DIR"
echo "  Cohort  : $COHORT"
echo "  Seeds   : ${SEEDS:-original}"
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
link_item "$REPO_DIR/config/bifo_mapping_ddkg.yaml" "$WORK_DIR/bifo_mapping_ddkg.yaml"

# Seed files
SEED_SUFFIX="${SEEDS:+_${SEEDS}}"
SEEDS_FILE="$REPO_DIR/data/cohorts/${COHORT}/kf_${COHORT}_seeds${SEED_SUFFIX}.txt"
link_item "$SEEDS_FILE" "$WORK_DIR/kf_${COHORT}_seeds${SEED_SUFFIX}.txt"
if [ -n "$SEEDS" ]; then
    link_item "$SEEDS_FILE" "$WORK_DIR/kf_${COHORT}_seeds.txt"
fi

echo ""
echo "Workspace ready: $WORK_DIR"
echo ""
echo "Next steps:"
echo "  cd $WORK_DIR"
echo "  bash $SCRIPT_DIR/run_full_pipeline.sh $COHORT ${SEEDS} neo4j PASSWORD bolt://localhost:7687"
