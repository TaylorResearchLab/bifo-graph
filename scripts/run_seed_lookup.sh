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

echo "============================================================"
echo "Seed CUI lookup for cohort: $COHORT"
echo "============================================================"

python3 "$SCRIPT_DIR/seed_cui_lookup.py" \
    --gene-list  "kf_${COHORT}_seeds${SEEDS:+_${SEEDS}}.txt" \
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
