#!/usr/bin/env bash
set -euo pipefail
COHORT="${1:-chd}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
python3 "$REPO_DIR/pipeline/build_cilia_reference.py" \
    --scores "kf_${COHORT}_results/pathway_scores_standard.csv" \
    --out    "kf_${COHORT}_cilia_reference.txt" \
    --verbose
echo "Next step: bash run_baseline.sh $COHORT cilia"
