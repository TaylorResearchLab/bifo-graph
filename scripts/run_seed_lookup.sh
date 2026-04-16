#!/usr/bin/env bash
# =============================================================================
# run_seed_lookup.sh
# Stage 2.4: Map variant-derived gene symbols to UMLS CUIs.
#
# Usage:
#   bash run_seed_lookup.sh [COHORT]
#   COHORT defaults to "chd"
#
# Inputs:
#   kf_{cohort}_seeds.txt         (gene symbols + carrier counts)
#   kf_{cohort}_nodes_clean.csv   (graph nodes with HGNC SAB)
#
# Outputs:
#   kf_{cohort}_seed_cuis.txt     (UMLS CUIs for matched seed genes)
# =============================================================================

set -euo pipefail

COHORT="${1:-chd}"
SEEDS="${2:-maf001}"   # maf001 (default/recommended) | maf01 | "" (original 56-gene list)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"
DATA_DIR="$REPO_DIR/data"
CONFIG_DIR="$REPO_DIR/config"

echo "============================================================"
echo "Seed CUI lookup for cohort: $COHORT"
echo "============================================================"

SEEDS_FILE="kf_${COHORT}_seeds${SEEDS:+_${SEEDS}}.txt"
if [ ! -f "$SEEDS_FILE" ]; then
    REPO_SEEDS="$DATA_DIR/cohorts/${COHORT}/$SEEDS_FILE"
    if [ -f "$REPO_SEEDS" ]; then
        echo "  (using repo seed file: $REPO_SEEDS)"
        SEEDS_FILE="$REPO_SEEDS"
    else
        echo "ERROR: Seed file not found: $SEEDS_FILE"
        echo "Run setup_workspace.sh first."
        exit 1
    fi
fi

python3 "$PIPELINE_DIR/seed_cui_lookup.py" \
    --gene-list  "$SEEDS_FILE" \
    --nodes-csv  "kf_${COHORT}_nodes_clean.csv" \
    --out        "kf_${COHORT}_seed_cuis.txt" \
    --verbose

N=$(grep -vc '^#' "kf_${COHORT}_seed_cuis.txt" 2>/dev/null || echo 0)
echo ""
echo "Seed lookup complete: $N CUIs written."
if [ "$N" -lt 10 ]; then
    echo "WARNING: Fewer than 10 seeds found. Check nodes_clean.csv has HGNC SAB rows."
fi
echo "Next step: bash run_conditioning.sh $COHORT"
