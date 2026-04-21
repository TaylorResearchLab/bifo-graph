#!/usr/bin/env bash
# =============================================================================
# merge_files.sh
# Stage 2.2: Merge cleaned CSV files using pandas (header-safe).
#
# Usage:
#   bash merge_files.sh [COHORT]
#   COHORT defaults to "chd"
#
# Inputs:
#   kf_{cohort}_nodes_clean.csv
#   kf_{cohort}_pathway_member_nodes_clean.csv
#   kf_{cohort}_edges_raw_clean.csv
#   kf_{cohort}_pathway_membership_edges_clean.csv
#
# Outputs:
#   kf_{cohort}_nodes_extended.csv   (nodes_clean + pathway_member_nodes)
#   kf_{cohort}_edges_all.csv        (edges_raw_clean + pathway_membership)
#
# NOTE: Uses pandas pd.concat — never cat. The cat command does not skip
# duplicate CSV headers and produces malformed files.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PIPELINE_DIR="$REPO_DIR/pipeline"

COHORT="${1:-chd}"

echo "============================================================"
echo "Merging CSV files for cohort: $COHORT"
echo "============================================================"

python3 - <<PYEOF
import pandas as pd

cohort = "${COHORT}"

def merge(files, out):
    dfs = [pd.read_csv(f, dtype=str) for f in files]
    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(out, index=False)
    print(f"  {out}: {len(merged):,} rows, columns={list(merged.columns)}")

merge(
    [f"kf_{cohort}_nodes_clean.csv",
     f"kf_{cohort}_pathway_member_nodes_clean.csv"],
    f"kf_{cohort}_nodes_extended.csv"
)

merge(
    [f"kf_{cohort}_edges_raw_clean.csv",
     f"kf_{cohort}_pathway_membership_edges_clean.csv"],
    f"kf_{cohort}_edges_all.csv"
)

print("")
print("Merge complete.")
PYEOF

echo "Next step: bash run_seed_lookup.sh $COHORT"
